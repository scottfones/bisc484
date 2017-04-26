# FLASK IMPORTS
from flask import render_template, flash, redirect, request, make_response
from app import app
from .forms import ProteinInputForm, ABIInputForm
from werkzeug import secure_filename


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
from Bio import Entrez, AlignIO, SeqIO, Phylo, SubsMat, Alphabet


# Define File Paths
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
ABS_TMP = os.path.join(APP_ROOT, 'static/tmp/')


# Config for File Uploads
UPLOAD_FOLDER = 'static/tmp/' 
ALLOWED_EXTENSIONS = set(['txt', 'fasta', 'faa', 'abi'])


# Check if uploaded file is of allowed type
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/dnaInput', methods=['GET', 'POST'])
def dnaInput():

    abi_form = ABIInputForm()

    if abi_form.validate_on_submit():

        cookieExpireDate = datetime.datetime.now()
        cookieExpireDate = cookieExpireDate + datetime.timedelta(days=90)

        abiRedirect = make_response(redirect('/abiResults'))

        userID = request.cookies.get('uuid')

        if userID is None:
            abiRedirect.set_cookie('uuid', str(uuid.uuid4()), expires=cookieExpireDate)

        abi_filename = ABS_TMP + userID + '.abi'
        abi_form.abi_file.data.save(abi_filename)

        return abiRedirect 

    return render_template('dnaInput.html', 
                           title = 'DNA Input',
                           abi_form = abi_form)




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


@app.route('/abiResults')
def abiResults():

    # Read User UUID
    userID = request.cookies.get('uuid')

    # Check History 
    if userID is None:
       flash(u'ABI file not found. Please resubmit.', 'error')
       return redirect('/dnaInput')

    abi_file = ABS_TMP + userID + '.abi'
    record = SeqIO.read(abi_file, 'abi')

    channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']

    trace = collections.defaultdict(list)

    for c in channels:
        trace[c] = record.annotations['abif_raw'][c]

    seq_length = len(trace['DATA9'])
    x = np.linspace(1,seq_length,seq_length)

    res_locations = record.annotations['abif_raw']['PLOC2']
    res_labels_string = record.annotations['abif_raw']['PBAS2']
    res_labels = []

    for char in res_labels_string:
        res_labels.append(char)

    label_y = []

    for location in res_locations:
        label_y.append(max(trace['DATA9'][location], trace['DATA10'][location], trace['DATA11'][location], trace['DATA12'][location]) + 50)

    quality_scores =  record.letter_annotations['phred_quality']

    dye_label = []

    for i in range(1,5):
        dye = 'DyeW' + str(i)
        wavelength = record.annotations['abif_raw'][dye]

        if wavelength == 540:
            dye_label.append('Wavelength 540 (G)')
        elif wavelength == 568:
            dye_label.append('Wavelength 568 (A)')
        elif wavelength == 595:
            dye_label.append('Wavelength 595 (T)')
        elif wavelength == 615:
            dye_label.append('Wavelength 615 (C)') 

    trim_sequence_record = SeqIO.AbiIO._abi_trim(record)

    raw_sequence = res_labels_string
    trim_sequence = trim_sequence_record.seq

    trace_1 = go.Scatter(
        x = x,
        y = trace['DATA9'],
        fill = 'tozeroy',
        line=dict(
            color='#440154',
        ),
        opacity = 0.7,
        name = dye_label[0],
    )
    trace_2 = go.Scatter(
        x = x,
        y = trace['DATA10'],
        fill = 'tozeroy',
        line=dict(
            color='#355F8D',
        ),
        opacity = 0.7,
        name = dye_label[1],
    )
    trace_3 = go.Scatter(
        x = x,
        y = trace['DATA11'],
        fill = 'tozeroy',
        line=dict(
            color='#44BE70',
        ),
        opacity = 0.7,
        name = dye_label[2],
    )
    trace_4 = go.Scatter(
        x = x,
        y = trace['DATA12'],
        fill = 'tozeroy',
        line=dict(
            color='#FDE725',
        ),
        opacity = 0.8,
        name = dye_label[3],
    )
    trace_seq = go.Scatter(
        x = res_locations,
        y = label_y,
        mode = 'text',
        name = 'Basecaller Sequence',
        marker = dict(
            size = 3,
        ),
        text = res_labels,
    )
    trace_quality = go.Scatter(
        x = res_locations,
        y = quality_scores,
        name = 'Quality Scores',
        line = dict(
            width = 5,
            color = '#dd4444',
        ),
        yaxis = 'y2',
    )
    
    data = [trace_1,trace_2,trace_3,trace_4,trace_seq,trace_quality]

    layout = go.Layout (

        title = 'ABI File Raw Plots',
        hovermode= 'closest',
        showlegend = True,
        width = 800,
        height = 700,
        autosize = False,

        font = dict(
            size = 22,
        ),

        xaxis = dict(
            title = 'sequence data',
            range = [0,200],
            showgrid = False,
            titlefont=dict(
                size=20,
            ),
            tickfont=dict(
                size=14,
            ),
            zerolinewidth=0,
            ticks = 'outside',
        ),

        yaxis = dict(
            title = 'signal strength',
            range = [0,2800],
            showgrid = False,
            titlefont=dict(
                size=20,
            ),
            tickfont=dict(
                size=14,
            ),
            zerolinewidth=0,
            ticks = 'outside',
        ),

        yaxis2=dict(
            title='quality score',
            range = [-50,50],
            titlefont=dict(
                color='#aa1111',
                size = 20,
            ),
            tickfont=dict(
               color='#aa1111',
               size = 14,
            ),
            zerolinewidth=0,
            overlaying='y',
            side='right'
        ),

        legend = dict(
            x = 0,
            y = 1,
            bordercolor = '#333333',
            borderwidth = 1,
            font = dict(
                size = 16
            ),
            xanchor = 'left',
            yanchor = 'top',
        ),
    )

    fig = go.Figure(data=data, layout=layout)

    abi_plot_div = plot(fig, output_type='div')
    
    return render_template('abiResults.html',
                            title = 'ABI Results',
                            raw_sequence = raw_sequence,
                            trim_sequence = trim_sequence,
                            abi_plot_div = abi_plot_div)


@app.route('/proteinResults')
def proteinResults():
    
    # Read User UUID
    userID = request.cookies.get('uuid')

    # Check History State
    if userID is None:
       flash(u'Records not found. Please resubmit.', 'error') 
       return redirect('/proteinInput')

    # Initialize sequence list
    seqFile = ABS_TMP + userID + '.faa'
    sequenceList = []

    accessionList = request.cookies.get('accessionValues')
    
    if accessionList ==  'file':
        with open(seqFile, 'r') as fileInput:
            for record in SeqIO.parse(fileInput, 'fasta'):
                sequenceList.append(record)

        if not os.path.isfile(seqFile):
            f = open(seqFile, 'w')

            for record in sequenceList:
                f.write('> ' + record.name + '\n' + str(record.seq) + '\n')

            f.close()

    else:
        accessionList = accessionList.split('+')

        if not os.path.isfile(seqFile):
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
            f = open(seqFile, 'w')

            for record in sequenceList:
                f.write('> ' + record.name + '\n' + str(record.seq) + '\n')

            f.close() 
        else:
            with open(seqFile, 'r') as fileInput:
                for record in SeqIO.parse(fileInput, 'fasta'):
                    sequenceList.append(record)

       


    # Create Alignment and Tree File
    alignFile = ABS_TMP + userID + '.afa'
    treeFile = ABS_TMP + userID + '.dnd'

    if not os.path.isfile(alignFile):
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
    texFile = ABS_TMP + userID + '_identity.tex'
    pdfFile = ABS_TMP + userID + '_identity.pdf'
    pngFile = ABS_TMP + userID + '_identity.png'

    if not os.path.isfile(pngFile):
        with open(texFile, 'w') as texfile:
        
            texfile.write('\\documentclass[preview]{standalone}\n\\usepackage{texshade}\n\\begin{document}\n\\begin{texshade}{%s}\n\\shadingmode[allmatchspecial]{identical}\n\\nomatchresidues{Gray70}{White}{upper}{bf}\n\\conservedresidues{Black}{LightCyan}{upper}{bf}\n\\allmatchresidues{White}{Red}{upper}{bf}\n\\showcaption[bottom]{\\textbf{Protein MSA with Identity Highlighting}}\n\\hideconsensus\n\\namesfootnotesize\n\\residuesfootnotesize\n\\numberingtiny\n\\legendfootnotesize\n\\showlegend\n\\end{texshade}\n\\end{document}\n' % alignFile)

        texfile.close()

        os.system('pdflatex -output-directory=%s %s' % (ABS_TMP,texFile))
        os.system('pdftoppm -r 250 -png %s > %s' % (pdfFile,pngFile))

    pngIdentityFile = 'static/tmp/' + userID + '_identity.png'


    #Generate Latex Chemical Alignment
    texFile = ABS_TMP + userID + '_chemical.tex'
    pdfFile = ABS_TMP + userID + '_chemical.pdf'
    pngFile = ABS_TMP + userID + '_chemical.png'

    if not os.path.isfile(pngFile):
        with open(texFile, 'w') as texfile:

            texfile.write('\\documentclass[preview]{standalone}\n\\usepackage{texshade}\n\\begin{document}\n\\begin{texshade}{%s}\n\\shadingmode[chemical]{functional}\n\\showcaption[bottom]{\\textbf{Protein MSA with Chemical Highlighting}}\n\\hideconsensus\n\\namesfootnotesize\n\\residuesfootnotesize\n\\numberingtiny\n\\legendfootnotesize\n\\showlegend\n\\end{texshade}\n\\end{document}\n' % alignFile)

        texfile.close()

        os.system('pdflatex -output-directory=%s %s' % (ABS_TMP,texFile))
        os.system('pdftoppm -r 250 -png %s > %s' % (pdfFile,pngFile))

    pngChemicalFile = 'static/tmp/' + userID + '_chemical.png'


    #Generate Latex Chemical Alignment
    texFile = ABS_TMP + userID + '_structural.tex'
    pdfFile = ABS_TMP + userID + '_structural.pdf'
    pngFile = ABS_TMP + userID + '_structural.png'

    if not os.path.isfile(pngFile):
        with open(texFile, 'w') as texfile:

            texfile.write('\\documentclass[preview]{standalone}\n\\usepackage{texshade}\n\\begin{document}\n\\begin{texshade}{%s}\n\\shadingmode[structure]{functional}\n\\showcaption[bottom]{\\textbf{Protein MSA with Structural Highlighting}}\n\\hideconsensus\n\\namesfootnotesize\n\\residuesfootnotesize\n\\numberingtiny\n\\legendfootnotesize\n\\showlegend\n\\end{texshade}\n\\end{document}\n' % alignFile)

        texfile.close()

        os.system('pdflatex -output-directory=%s %s' % (ABS_TMP,texFile))
        os.system('pdftoppm -r 250 -png %s > %s' % (pdfFile,pngFile))

    pngStructuralFile = 'static/tmp/' + userID + '_structural.png'


    # Create ASCII Dendrogram
    tree = Phylo.read(treeFile, "newick")
   
    # re-root tree 
    leghemeClade = tree.find_clades({'name':'.*leghemoglobin.*'}) 

    tree.root_with_outgroup(leghemeClade)

    leghemeClade = tree.find_clades({'name':'.*outgroup.*'}) 

    tree.root_with_outgroup(leghemeClade)

    # find and categorize term nodes
    terminalClades = tree.get_terminals()
    alphaClades = []
    betaClades = []
    cytoClades = []
    myoClades = []
    neuroClades = []

    for clade in terminalClades:
        if 'alpha' in clade.name.lower():
            alphaClades.append(clade.name) 
        if 'theta' in clade.name.lower():
            alphaClades.append(clade.name) 
        if 'mu' in clade.name.lower():
            alphaClades.append(clade.name) 
        if 'zeta' in clade.name.lower():
            alphaClades.append(clade.name) 
        if 'beta' in clade.name.lower():
            betaClades.append(clade.name) 
        if 'gamma' in clade.name.lower():
            betaClades.append(clade.name) 
        if 'delta' in clade.name.lower():
            betaClades.append(clade.name) 
        if 'epsilon' in clade.name.lower():
            betaClades.append(clade.name) 
        if 'cyto' in clade.name.lower():
            cytoClades.append(clade.name)
        if 'myo' in clade.name.lower():
            myoClades.append(clade.name)
        if 'neuro' in clade.name.lower():
            neuroClades.append(clade.name)

    for clade in alphaClades:
        mrca = tree.common_ancestor({'name':alphaClades[0]}, {'name':clade})
        mrca.color = 'red'

    for clade in betaClades:
        mrca = tree.common_ancestor({'name':betaClades[0]}, {'name':clade})
        mrca.color = 'blue'

    for clade in cytoClades:
        mrca = tree.common_ancestor({'name':cytoClades[0]}, {'name':clade})
        mrca.color = '#ff9900'

    for clade in myoClades:
        mrca = tree.common_ancestor({'name':myoClades[0]}, {'name':clade})
        mrca.color = 'brown'

    for clade in neuroClades:
        mrca = tree.common_ancestor({'name':neuroClades[0]}, {'name':clade})
        mrca.color = 'green'


    # Create Graphic Dendrogram
    graphicTreeFile = ABS_TMP + userID + '_graphicTree.png'

    if not os.path.isfile(graphicTreeFile):

        tree = tree.as_phyloxml()
        tree.root.color='gray'

        matplotlib.rc('font', size=16)
        matplotlib.rc('lines', linewidth=4.0)
        pylab.rcParams['figure.figsize'] = 12.5, 10
        pylab.rcParams['figure.autolayout'] = True
        pylab.rcParams['savefig.bbox'] = 'tight'
        Phylo.draw(tree, do_show=False, axes=None)
        pylab.axis('off')
        
        pylab.savefig(graphicTreeFile, dpi=150)

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
                    font = dict(color='#E0E0E0'if val < 5 else '#222222', size=12),
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

    scoreTable.append('<table class="table table-sm table-responsive"><tr><td></td>')

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

 
    # Generate Substitution Matrix by Residue Type
    colorMat = np.array([[.1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.0,.2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.9,.2,.2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.0,.2,.2,.2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.9,.0,.9,.0,.3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.0,.9,.0,.9,.3,.3,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.9,.0,.9,.0,.9,.0,.4,1,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.0,.9,.0,.9,.0,.9,.4,.4,1,1,1,1,1,1,1,1,1,1,1,1],
                         [.9,.0,.9,.0,.9,.0,.4,.4,.4,1,1,1,1,1,1,1,1,1,1,1],
                         [.0,.9,.0,.9,.0,.9,.4,.4,.4,.4,1,1,1,1,1,1,1,1,1,1],
                         [.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.5,1,1,1,1,1,1,1,1,1],
                         [.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.6,1,1,1,1,1,1,1,1],
                         [.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.6,.6,1,1,1,1,1,1,1],
                         [.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.6,.6,.6,1,1,1,1,1,1],
                         [.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.6,.6,.6,.6,1,1,1,1,1],
                         [.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.6,.6,.6,.6,.6,1,1,1,1],
                         [.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.7,1,1,1],
                         [.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.7,.7,1,1],
                         [.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.7,.7,.7,1],
                         [.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.9,.0,.8]])

    resTypeKeys = ['C', 'R', 'K', 'H', 'E', 'D', 'Q', 'N', 'T', 'S', 
                   'G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'P'] 
            
    resTypeAnn = []

    for x_i, xPro in enumerate(resTypeKeys):
        rowList = []

        for y_i in range(0,(x_i+1)):
            yPro = resTypeKeys[y_i]

            transX = proteinKeys.index(xPro) 
            transY = proteinKeys.index(yPro) + 1

            rowList.append(int(scoreTxtMat[transX][transY]))

        resTypeAnn.append(rowList)

    hover = []

    for x_i, xPro in enumerate(resTypeKeys):
        rowList = []

        for y_i in range(0,(x_i+1)):
            yPro = resTypeKeys[y_i]

            xProType = ''
            yProType = ''

            if xPro == 'C':
                xProType = 'Hydrophilic, Special Case'
            elif xPro in ['R', 'K', 'H']:
                xProType = 'Hydrophilic, Basic'
            elif xPro in ['E', 'D']:
                xProType = 'Hydrophilic, Acidic'
            elif xPro in ['Q','N','T','S']:
                xProType = 'Hydophilic, Polar'
            elif xPro == 'G':
                xProType = 'No Side Chain'
            elif xPro in ['A','V','I','L','M']:
                xProType = 'Hydrophobic, Aliphatic'
            elif xPro in ['F','Y','W']:
                xProType = 'Hydrophobic, Aromatic'
            elif xPro == 'P':
                xProType = 'Hydrophobic, Special Case'

            if yPro == 'C':
                yProType = 'Hydrophilic, Special Case'
            elif yPro in ['R', 'K', 'H']:
                yProType = 'Hydrophilic, Basic'
            elif yPro in ['E', 'D']:
                yProType = 'Hydrophilic, Acidic'
            elif yPro in ['Q','N','T','S']:
                yProType = 'Hydophilic, Polar'
            elif yPro == 'G':
                yProType = 'No Side Chain'
            elif yPro in ['A','V','I','L','M']:
                yProType = 'Hydrophobic, Aliphatic'
            elif yPro in ['F','Y','W']:
                yProType = 'Hydrophobic, Aromatic'
            elif yPro == 'P':
                yProType = 'Hydrophobic, Special Case'
            
            rowList.append(xPro + ', ' + yPro + ': ' + str(resTypeAnn[x_i][y_i]) + '<br>' + xPro + ': ' + xProType + '<br>' + yPro + ': ' + yProType)
         
        hover.append(rowList)


    annotations = []

    for n, row in enumerate(resTypeAnn):
        for m, val in enumerate(row):
            var = resTypeAnn[n][m]
            
            resCat = colorMat[n][m]

            if resCat == .0 or resCat == .9:
                preColor = '#101010'
            if resCat == .1:
                preColor = '#FAFAFA'
            elif resCat > .1 and resCat < .5:
                preColor = '#E0E0E0'
            elif resCat == .5:
                preColor = '#D0D0D0'
            elif resCat > .6 and resCat < .8:
                preColor = '#303030'
            elif resCat == .8:
                preColor = '#000000'


            annotations.append(
                dict(
                    text = '%i' % val,
                    x = resTypeKeys[m],
                    y = resTypeKeys[n],
                    xref = 'x1', yref= 'y1',
                    font = dict(
                       color = preColor, 
                       size=12),
                    showarrow = False) 
                )

    colorscale =    [[0.0, '#F9F9F9'], [.1, '#440154'], [.2, '#414387'], [.3, '#355F8D'], [.4, '#21918C'],
                     [.5, '#24A884'], [.6, '#44BE70'], [.7, '#7BD151'], [.8, '#FDE725'],  [.9, '#F5F5F5'],
                     [1.0, '#FFFFFF']]

    trace = go.Heatmap(x=resTypeKeys, y=resTypeKeys, z=colorMat, text=hover, hoverinfo='text', colorscale=colorscale, showscale=False)

    fig = go.Figure(data=[trace])
    fig['layout'].update(
        #title = 'Substitution Matrix Heatmap<br>Residues By Type',
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
            t = 0)
    )

    resTypeDiv = plot(fig, output_type='div')

    

    return render_template('proteinResults.html',
                            title='Protein Results',
                            sequenceList = sequenceList,
                            alignmentClustal = alignmentClustal,
                            alignmentFASTA = alignmentFASTA,
                            pngIdentityFile = pngIdentityFile,
                            pngChemicalFile = pngChemicalFile,
                            pngStructuralFile = pngStructuralFile,
                            graphicTree = graphicTreeFile,
                            singleLtrFrqDiv = singleLtrFrqDiv,
                            pairFrqHtMpDiv = pairFrqHtMpDiv,
                            scoreMatDiv = scoreMatDiv,
                            resTypeDiv = resTypeDiv,
                            scoreTable = scoreTable)




@app.route('/')
@app.route('/index')
def index():

    return render_template('index.html',
                           title = 'Home')
