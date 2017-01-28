# FLASK IMPORTS
from flask import render_template, flash, redirect, request, make_response
from app import app
from .forms import LoginForm, ProteinInputForm
from werkzeug.utils import secure_filename


# 'PYTHON' IMPORTS
import os
import glob
import math
import uuid
import datetime
import matplotlib
matplotlib.use('Agg')
import pylab
import numpy as np
import collections
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from decimal import Decimal
from plotly.offline import download_plotlyjs, plot
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Entrez, AlignIO, SeqIO, Phylo, SubsMat, Alphabet


# Define File Paths
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
ABS_TMP = os.path.join(APP_ROOT, 'static/tmp/')


# Config for File Uploads
UPLOAD_FOLDER = 'static/tmp/' 
ALLOWED_EXTENSIONS = set(['txt', 'fasta', 'faa'])


# Check if uploaded file is of allowed type
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/proteinInput', methods=['GET', 'POST'])
def proteinInput():

    form1 = ProteinInputForm()

    if form1.validate_on_submit():
    
        formList = form1.accessionInput.data.split( )

        cookieExpireDate = datetime.datetime.now()
        cookieExpireDate = cookieExpireDate + datetime.timedelta(days=90)

        accessionRedirect = make_response(redirect('/proteinResults'))
        accessionRedirect.set_cookie('accessionValues', '+'.join(formList), expires=cookieExpireDate)
       
        userID = request.cookies.get('uuid')

        if userID is None:
            accessionRedirect.set_cookie('uuid', str(uuid.uuid4()), expires=cookieExpireDate)
        else:
            for oldFile in glob.glob(ABS_TMP + userID + '*'):
                os.remove(oldFile)
            
            accessionRedirect.set_cookie('uuid', str(uuid.uuid4()), expires=cookieExpireDate)


        return accessionRedirect

    if request.method == 'POST':
        if 'file' not in request.files:
            flash('No file found')
            return redirect(request.url)

        file = request.files['file']

        if file.filename == '':
            flash('No file selected')
            return(redirect(request.url))

        if file and allowed_file(file.filename):

            cookieExpireDate = datetime.datetime.now()
            cookieExpireDate = cookieExpireDate + datetime.timedelta(days=90)

            accessionRedirect = make_response(redirect('/proteinResults'))
            accessionRedirect.set_cookie('accessionValues', 'file', expires=cookieExpireDate)

            userID = request.cookies.get('uuid')

            if userID is None:
                userID = str(uuid.uuid4())
                accessionRedirect.set_cookie('uuid', userID, expires=cookieExpireDate)
            else:
                for oldFile in glob.glob(ABS_TMP + userID + '*'):
                    os.remove(oldFile)

                userID = str(uuid.uuid4())
                accessionRedirect.set_cookie('uuid', userID, expires=cookieExpireDate)

            #filename = secure_filename(file.filename)

            file.save(ABS_TMP + userID + '.faa')

            
     
            return accessionRedirect



    if os.path.isfile(APP_ROOT + '/static/protein/accessionValues.txt'):

        inFile = open(APP_ROOT + '/static/protein/accessionValues.txt', 'r')
        accessionValues = []

        for inputLine in inFile:
            accessionValues.append(inputLine.strip()) 

        inFile.close()

        accessionList = "\n".join(accessionValues)

        form1 = ProteinInputForm(accessionInput = accessionList)

    return render_template('proteinInput.html',
                            title = 'Protein Input',
                            form1 = form1)



@app.route('/proteinResults')
def proteinResults():
    
    # Read User UUID
    userID = request.cookies.get('uuid')
    seqFile = ABS_TMP + userID + '.faa'
    sequenceList = []

    accessionList = request.cookies.get('accessionValues')
    
    if accessionList ==  'file':
        with open(seqFile, 'r') as fileInput:
            for record in SeqIO.parse(fileInput, 'fasta'):
                sequenceList.append(record)

        f = open(seqFile, 'w')

        for record in sequenceList:
            f.write('> ' + record.name + '\n')
            f.write(str(record.seq) + '\n')

        f.close()

    else:
        accessionList = accessionList.split('+')

        # Accession Look Up
        Entrez.email = "sfones@udel.edu"
        sequenceList = []

        for accessionValue in accessionList:
            try:
                record = SeqIO.read(Entrez.efetch(db="protein", id=accessionValue, rettype="fasta", retmode="text"), "fasta")
            except:
                flash(u'Accession Value Error', 'error')
                return redirect('/proteinInput')

            tmpDesc = record.description.split( )
            record.name = "_".join(tmpDesc[1:-2])
            record.description = "_".join(tmpDesc)

            sequenceList.append(record)


        # Construct Sequence List 
        alignmentList = []

        for record in sequenceList:

            alignmentList.append('>' + str(record.name))
            alignmentList.append(str(record.seq))

        alignmentInput = "\n".join(alignmentList)


        # Create Sequence File
        seqFile = ABS_TMP + userID + '.faa'

        f = open(seqFile, 'w')

        for record in sequenceList:

            f.write('> ' + record.name + '\n')
            f.write(str(record.seq) + '\n')

        f.close() 


    # Create Alignment and Tree File
    alignFile = ABS_TMP + userID + '.afa'
    treeFile = ABS_TMP + userID + '.dnd'

    clustalw_cline = ClustalwCommandline("clustalw2", infile=seqFile, outfile=alignFile, newtree=treeFile, outorder="aligned", align=True)
    
    try:
        clustalw_cline()
    except:
        flash(u'Alignment Error: Check formatting (Spaces in names? All names unique?)','error')
        return(redirect('/proteinInput'))
        
    
    # Generate and Format Alignments
    alphabet = Alphabet.Gapped(IUPAC.protein)
    alignment = AlignIO.read(alignFile, "clustal", len(sequenceList), alphabet)

    alignmentClustal = alignment.format('clustal')
    alignmentClustal = alignmentClustal.split('\n')

    alignmentFASTA = alignment.format('fasta')
    alignmentFASTA = alignmentFASTA.split('\n')

    
    # Generate Latex Identity Alignment
    texIdentityFile = ABS_TMP + userID + '_identity.tex'

    with open(texIdentityFile, 'w') as texfile:
       
        texfile.write('\\documentclass[preview]{standalone}\n')
        texfile.write('\\usepackage{texshade}\n')
        texfile.write('\\usepackage{inconsolata}\n')
        texfile.write('\\begin{document}\n')
        texfile.write('\\begin{texshade}{%s}\n' % alignFile)
        texfile.write('\\shadingmode[allmatchspecial]{identical}\n')
        texfile.write('\\showcaption[bottom]{\\textbf{Protein MSA with Identity Highlighting}}\n')
        texfile.write('\\label{fig:blast_tc66374}\n')
        texfile.write('\\hideconsensus\n')
        texfile.write('\\namesfootnotesize\n')
        texfile.write('\\residuesfootnotesize\n')
        texfile.write('\\numberingscriptsize\n')
        texfile.write('\\showlegend\n')
        texfile.write('\\end{texshade}\n')
        texfile.write('\\end{document}\n')

    texfile.close()

    os.system('pdflatex -output-directory=%s %s' % (ABS_TMP,texIdentityFile))
    
    texIdentityPDF = 'static/tmp/' + userID + '_identity.pdf'


    # Generate Latex Chemical Similarity Alignment
    texChemicalFile= ABS_TMP + userID + '_chemical.tex'

    with open(texChemicalFile, 'w') as texfile:
       
        texfile.write('\\documentclass[preview]{standalone}\n')
        texfile.write('\\usepackage{texshade}\n')
        texfile.write('\\usepackage{inconsolata}\n')
        texfile.write('\\begin{document}\n')
        texfile.write('\\begin{texshade}{%s}\n' % alignFile)
        texfile.write('\\shadingmode[chemical]{functional}\n')
        texfile.write('\\showcaption[bottom]{\\textbf{Protein MSA with Chemical Similarity Highlighting}}\n')
        texfile.write('\\hideconsensus\n')
        texfile.write('\\namesfootnotesize\n')
        texfile.write('\\residuesfootnotesize\n')
        texfile.write('\\numberingscriptsize\n')
        texfile.write('\\showlegend\n')
        texfile.write('\\end{texshade}\n')
        texfile.write('\\end{document}\n')

    texfile.close()

    os.system('pdflatex -output-directory=%s %s' % (ABS_TMP,texChemicalFile))
    
    texChemicalPDF = 'static/tmp/' + userID + '_chemical.pdf'    


    # Generate Latex Structural Similarity Alignment
    texStructuralFile = ABS_TMP + userID + '_structural.tex'

    with open(texStructuralFile, 'w') as texfile:
       
        texfile.write('\\documentclass[preview]{standalone}\n')
        texfile.write('\\usepackage{texshade}\n')
        texfile.write('\\usepackage{inconsolata}\n')
        texfile.write('\\begin{document}\n')
        texfile.write('\\begin{texshade}{%s}\n' % alignFile)
        texfile.write('\\shadingmode[structure]{functional}\n')
        texfile.write('\\showcaption[bottom]{\\textbf{Protein MSA with Structural Similarity Highlighting}}\n')
        texfile.write('\\hideconsensus\n')
        texfile.write('\\namesfootnotesize\n')
        texfile.write('\\residuesfootnotesize\n')
        texfile.write('\\numberingscriptsize\n')
        texfile.write('\\showlegend\n')
        texfile.write('\\end{texshade}\n')
        texfile.write('\\end{document}\n')

    texfile.close()

    os.system('pdflatex -output-directory=%s %s' % (ABS_TMP,texStructuralFile))
    
    texStructuralPDF = 'static/tmp/' + userID + '_structural.pdf'  


    # Create ASCII Dendrogram
    tree = Phylo.read(treeFile, "newick")
   
    # re-root tree 
    leghemeClade = tree.find_clades({'name':'.*leghemoglobin.*'}) 
    tree.root_with_outgroup(leghemeClade)
   
    alphaMRCA = tree.common_ancestor({'name':'.*zeta.*'}, {'name':'.*mu.*'})
    alphaMRCA.color = 'red'
    
    betaMRCA = tree.common_ancestor({'name':'.*epsilon.*'}, {'name':'.*delta.*'})
    betaMRCA.color = 'blue'

    cytoClades = tree.find_clades({'name':'.*cyto.*'})
    for clade in cytoClades:
        clade.color = '#ff9900'

    myoClades = tree.find_clades({'name':'.*myo.*'})
    for clade in myoClades:
        clade.color = 'brown'

    neuroClades = tree.find_clades({'name':'.*neuro.*'})
    for clade in neuroClades:
        clade.color = 'green'

    asciiTreeFile = ABS_TMP + userID + '_asciiTree.txt'
    asciiTree = []

    with open(asciiTreeFile, 'w') as asciiTF:
        Phylo.draw_ascii(tree, file = asciiTF)

    asciiTF.close()

    inFile = open(asciiTreeFile, 'r')

    for inLine in inFile:
        asciiTree.append(inLine)

    inFile.close()


    # Create Graphic Dendrogram
    graphicTreeFile = ABS_TMP + userID + '_graphicTree.png'

    tree = tree.as_phyloxml()
    tree.root.color='gray'

    matplotlib.rc('font', size=24)
    matplotlib.rc('lines', linewidth=4.0)
    pylab.rcParams['figure.figsize'] = 25, 20
    pylab.rcParams['figure.autolayout'] = True
    pylab.rcParams['savefig.bbox'] = 'tight'
    Phylo.draw(tree, do_show=False, axes=None)
    pylab.axis('off')
    
    pylab.savefig(graphicTreeFile, dpi=200)

    graphicTreeFile = 'static/tmp/' + userID + '_graphicTree.png'


    # Calculate Single Letter Frequencies
    proteinKeys = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
            
    totalLetters = 0
    resFrqs = dict.fromkeys(proteinKeys, 0)
  
    # get letter counts and total
    for record in alignment:
        for letter in record.seq:
            if (letter == '-'):
                continue
            else:
                resFrqs[letter] = resFrqs[letter] + 1
                totalLetters += 1
   
    # convert from counts to freqs
    # prime frqXY for bar graph
    frqX = []
    frqY = []

    for key in sorted(resFrqs):
        resFrqs[key] = resFrqs[key] / float(totalLetters)
        frqX.append(key)
        frqY.append(resFrqs[key])

    # Generate Single Frequency Table
    data = [go.Bar(
        x = frqX,
        y = frqY,

        marker = dict(
            color = 'rgb(58,79,134)',
            line = dict(
                color = 'rgb(8,48,107)',
                width = 1.5),
        ),
            #opacity = 0.7
    )]

    layout = go.Layout(
        title = 'Single Letter Frequencies',
        autosize = False,
        width = 700,
        height = 500,

        annotations = [
            dict(
                x = xi,
                y = yi,
                text = '%.1E' % Decimal(yi),
                font = dict(size=8),
                xanchor = 'center',
                yanchor = 'bottom',
                showarrow = False,
            ) for xi,yi in zip(frqX,frqY)
        ],

        xaxis = dict(
            title = 'Residue',
            #showgrid = True,
        ),

        yaxis = dict(
            title = 'Frequency',
            #showgrid = True,
        ),

        margin = dict(
            l = 55,
            r = 0,
            b = 45,
            t = 45
        )

    )

    fig = go.Figure(data=data, layout=layout)

    singleLtrFrqDiv = plot(fig, output_type='div') 


    # Calculate Pairwise Probabilities
    alignSummary = AlignInfo.SummaryInfo(alignment)
    
    replaceInfo = alignSummary.replacement_dictionary()

    accRepMat = SubsMat.SeqMat(replaceInfo)
    obsFrqMat = SubsMat._build_obs_freq_mat(accRepMat)
 
    obsFrq_z = []

    for x_i, xPro in enumerate(proteinKeys):
        rowList = []

        for y_i in range(0,(x_i+1)):
            rowList.append(obsFrqMat[(proteinKeys[y_i],xPro)])

        obsFrq_z.append(rowList)


    # Generate Pairwise Probability Heatmap
    annotations = []

    for n, row in enumerate(obsFrq_z):
        for m, val in enumerate(row):
            var = obsFrq_z[n][m]

            annotations.append(
                dict(
                    text = '%.1E' % Decimal(val) if val != 0 else '0',
                    x = proteinKeys[m],
                    y = proteinKeys[n],
                    xref = 'x1', yref= 'y1',
                    font = dict(color='#E0E0E0'if val < 0.035 else 'black', size=8),
                    showarrow = False)
            )
  
    colorscale = 'Viridis'
    
    trace = go.Heatmap(x=proteinKeys, y=proteinKeys, z=obsFrq_z, colorscale=colorscale, showscale=False)

    fig = go.Figure(data=[trace])
    fig['layout'].update(
        title = "Pairwise Probability Heatmap",
        annotations=annotations,
        width = 700,
        height = 700,
        autosize = False,

        xaxis = dict(
            range = [-0.5,19.5],
            autorange =  True,
            type = 'category',
            ticks='', side='bottom',
            showgrid=False,
        ),
        yaxis = dict(
            range = [19.5,-0.5],
            autorange = True,
            showgrid=False,
            ticks='', ticksuffix='  ',
            type = 'category',
        ),

        margin = dict(
            l = 25,
            r = 0,
            b = 45,
            t = 45
        )
    )

    pairFrqHtMpDiv = plot(fig, output_type='div')


    # Calculate Substitution Matrix
    lambdaValue = 0.347
    scoreMat_z = []

    for x_i, xPro in enumerate(proteinKeys):
        rowList = []

        for y_i in range(0,(x_i+1)):
            if ( obsFrqMat[(proteinKeys[y_i],xPro)] == 0 ):
                rowList.append(-4)    ########### MIN VALUE THRESHOLD
            else:
                newSubValue = round( (1/lambdaValue) * np.log( obsFrqMat[(proteinKeys[y_i],xPro)] / (resFrqs[xPro] * resFrqs[proteinKeys[y_i]]) ))

                if ( newSubValue < -4 ):    ############ MIN VALUE THRESHOLD
                   newSubValue = -4 

                rowList.append(newSubValue)


        scoreMat_z.append(rowList)
    
    # Generate Scoring Matrix Heatmap
    annotations = []

    for n, row in enumerate(scoreMat_z):
        for m, val in enumerate(row):
            var = scoreMat_z[n][m]

            annotations.append(
                dict(
                    text = '%i' % val,
                    x = proteinKeys[m],
                    y = proteinKeys[n],
                    xref = 'x1', yref= 'y1',
                    font = dict(color='#E0E0E0'if val < 5 else 'black', size=12),
                    showarrow = False)
            )

    colorscale = 'Viridis'

    trace = go.Heatmap(x=proteinKeys, y=proteinKeys, z=scoreMat_z, colorscale=colorscale, showscale=False)

    fig = go.Figure(data=[trace])
    fig['layout'].update(
        title = 'Substitution Matrix Heatmap',
        annotations=annotations,
        width = 700,
        height = 700,
        autosize = False,

        xaxis = dict(
            range = [-0.5,19.5],
            autorange =  True,
            type = 'category',
            ticks='', side='bottom',
            showgrid=False,
            #title = 'Residue 1'
        ),
        yaxis = dict(
            range = [19.5,-0.5],
            autorange = True,
            showgrid=False,
            ticks='', ticksuffix='  ',
            type = 'category',
            #title = 'Residue 2'
        ),

        margin = dict(
            l = 25,
            r = 0,
            b = 45,
            t = 45
        )
    )

    scoreMatDiv = plot(fig, output_type='div')


    # Translate Scoring Matrix to Text
    #scoreTxtMat = np.zeros((20,20))
    #scoreTxtMat = scoreTxtMat + np.array(scoreMat_z)

    for i,row in enumerate(scoreMat_z):
        scoreMat_z[i].extend(np.zeros(20-i))

    scoreTxtMat = np.array(scoreMat_z)

    sm_rows,sm_cols = np.triu_indices(len(scoreTxtMat),1)
    
    scoreTxtMat[sm_rows,sm_cols] = scoreTxtMat[sm_cols,sm_rows]

    scoreTxtMat = np.delete(scoreTxtMat,20,1)
    scoreTxtMat = scoreTxtMat.tolist()

    for i,row in enumerate(scoreTxtMat):
        for j,num in enumerate(row):
            scoreTxtMat[i][j] = str(int(num))
        scoreTxtMat[i].insert(0,proteinKeys[i])

    scoreTable = []

    scoreTable.append('<table class="table table-condensed table-responsive"><tr><td></td>')
    
    for letter in proteinKeys:
        scoreTable.append('<td>' + letter + '</td>')

    scoreTable.append('</tr>')

    for row in scoreTxtMat:
        scoreTable.append('<tr>')
        for element in row:
            scoreTable.append('<td>' + element + '</td>')
        scoreTable.append('</tr>')

    scoreTable.append('</table>')
    scoreTable = ' '.join(scoreTable)



    return render_template('proteinResults.html',
                            title='Protein Results',
                            sequenceList = sequenceList,
                            alignmentClustal = alignmentClustal,
                            alignmentFASTA = alignmentFASTA,
                            texIdentityPDF = texIdentityPDF,
                            texChemicalPDF = texChemicalPDF,
                            texStructuralPDF = texStructuralPDF,
                            asciiTree = asciiTree,
                            graphicTree = graphicTreeFile,
                            singleLtrFrqDiv = singleLtrFrqDiv,
                            pairFrqHtMpDiv = pairFrqHtMpDiv,
                            scoreMatDiv = scoreMatDiv,
                            scoreTable = scoreTable)



@app.route('/')
@app.route('/index')
def index():

    return render_template('index.html',
                           title='Home')
