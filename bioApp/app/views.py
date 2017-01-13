# FLASK IMPORTS
from flask import render_template, flash, redirect, request, make_response
from app import app
from .forms import LoginForm, ProteinInputForm

# 'PYTHON' IMPORTS
import os
import glob
import uuid
import datetime
import matplotlib
matplotlib.use('Agg')
import pylab
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from plotly.offline import download_plotlyjs, plot
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Entrez, AlignIO, SeqIO, Phylo, SubsMat, Alphabet


# Define File Paths
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
ABS_TMP = os.path.join(APP_ROOT, 'static/tmp/')


'''
@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    
    if form.validate_on_submit():
        flash('Login requested for OpenID="%s", remember_me=%s' %
              (form.openid.data, str(form.remember_me.data)))
        return redirect('/index')
   
    return render_template('login.html', 
                            title='Sign In',
                            form=form,
                            providers=app.config['OPENID_PROVIDERS'])
'''

@app.route('/proteinInput', methods=['GET', 'POST'])
def proteinInput():

    form = ProteinInputForm()

    if form.validate_on_submit():
    
        formList = form.accessionInput.data.split( )

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

    if os.path.isfile(APP_ROOT + '/static/protein/accessionValues.txt'):

        inFile = open(APP_ROOT + '/static/protein/accessionValues.txt', 'r')
        accessionValues = []

        for inputLine in inFile:
            accessionValues.append(inputLine.strip()) 

        inFile.close()

        accessionList = "\n".join(accessionValues)

        form = ProteinInputForm(accessionInput = accessionList)

    return render_template('proteinInput.html',
                            title='Protein Input',
                            form=form)


@app.route('/proteinResults')
def proteinResults():

    # Read and Format Accession Values from Cookie
    accessionList = request.cookies.get('accessionValues')
    accessionList = accessionList.split('+')

    # Read User UUID
    userID = request.cookies.get('uuid')

    # Accession Look Up
    Entrez.email = "sfones@udel.edu"
    sequenceList = []

    for accessionValue in accessionList:
    
        record = SeqIO.read(Entrez.efetch(db="protein", id=accessionValue, rettype="fasta", retmode="text"), "fasta")
        
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

#    clustalomega_cline = ClustalOmegaCommandline(infile=seqFile, outfile=alignFile, guidetree_out=treeFile, outfmt="clu", outputorder='tree-order', force=True)
#    clustalomega_cline()


    clustalw_cline = ClustalwCommandline("clustalw2", infile=seqFile, outfile=alignFile, newtree=treeFile, outorder="aligned", align=True)
    clustalw_cline()

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
    pylab.rcParams['figure.figsize'] = 15, 10
    pylab.rcParams['figure.autolayout'] = True
    pylab.rcParams['savefig.bbox'] = 'tight'
    Phylo.draw(tree, do_show=False, axes=None)
    pylab.axis('off')
    
    pylab.savefig(graphicTreeFile, dpi=200)

    graphicTreeFile = 'static/tmp/' + userID + '_graphicTree.png'


    # Calculate Alignment Statistics
    alignSummary = AlignInfo.SummaryInfo(alignment)

    replaceInfo = alignSummary.replacement_dictionary()

    proteinReplacementCounts_x = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                                  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    proteinReplacementCounts_y = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                                  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    replace_x = []
    replace_y = []
    replace_z = []
    
    for x_i, xPro in enumerate(proteinReplacementCounts_x):
        for y_i, yPro in enumerate(proteinReplacementCounts_y):
            replace_x.append(xPro)
            replace_y.append(yPro)
            replace_z.append(replaceInfo[(xPro,yPro)])

    # Generate Acceptable Count Matrix
    #accRepMat = SubsMat.SeqMat(replaceInfo)
    #accRepMat.print_full_mat()

    accRep_x = []
    accRep_y = []
    accRep_z = []

    for x_i, xPro in enumerate(proteinReplacementCounts_x):
        for y_i, yPro in enumerate(proteinReplacementCounts_y):
            accRep_x.append(xPro)
            accRep_y.append(yPro)
            
            if (xPro == yPro):
                accRep_z.append(replaceInfo[(xPro,yPro)])
            else:
                accRep_z.append(replaceInfo[(xPro,yPro)] + replaceInfo[(yPro,xPro)])

    # Generate Accepted Replacement Table
    accRepTableHTML = []

    accRepTableHTML.append('<table style="width:700px;table-layout:fixed" class="table table-condensed table-striped table-bordered">')
    accRepTableHTML.append('<thead><tr><th></th>')

    for residue in proteinReplacementCounts_x:
        accRepTableHTML.append('<th> %s </th>' % residue)

    accRepTableHTML.append('</tr></thead><tbody>')

    for residuex in proteinReplacementCounts_x:
        accRepTableHTML.append('<tr><th> %s </th>' % residuex)
        
        for residuey in proteinReplacementCounts_y:
            if (residuex == residuey):
                accRepTableHTML.append('<td> %g </td>' % replaceInfo[(residuey,residuex)])
            else:
                repSum = replaceInfo[(residuey,residuex)] + replaceInfo[(residuex,residuey)]
                accRepTableHTML.append('<td> %g </td>' % repSum)

        accRepTableHTML.append('</tr>')

    accRepTableHTML.append('</tbody></table>')

    accRepTable = ''.join(accRepTableHTML)


    # Generate 3D Accepted Replacement Plot
    trace1 = go.Scatter3d(
        x = accRep_x,
        y = accRep_y,
        z = accRep_z,
       
        mode = 'markers',

        marker = dict(
            size = 6,
            color = accRep_z,
            symbol = 'circle-dot',
            colorscale = 'Viridis',
            opacity = 0.8
        )
    )
    
    data = [trace1]
    layout = go.Layout(
        autosize = False,
        width = 700,
        height = 500,

        title = 'Accepted Replacement Counts',

        scene = dict(
            xaxis = dict(
                title = 'Residue 1',
                showgrid = True,
                dtick = 1
            ),

            yaxis = dict(
                title = 'Residue 2',
                showgrid = True,
                dtick = 1
            ),

            zaxis = dict(
                title = 'Count',
                showgrid = True,
                dtick = 100
            ),
        ),

        margin = dict(
            l = 0,
            r = 0,
            b = 0,
            t = 30
        )
    )

    fig = go.Figure(data=data, layout=layout)
    
    accRepName = ABS_TMP + userID + '_accRepPlot'
    
    accRepPlotDiv = plot(fig, output_type='div')
    
    # Calculate Log Odds Matrix
    baseE = 2.71828183 # e approx
    scalingFactor = 2.88184438 # 1 / lambda approx

    my_arm = SubsMat.SeqMat(replaceInfo)
    my_lom = SubsMat.make_log_odds_matrix(my_arm)



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
                            #replacementTable = replacementTable,
                            #replacePlotDiv = replacePlotDiv,
                            accRepTable = accRepTable,
                            accRepPlotDiv = accRepPlotDiv)



@app.route('/')
@app.route('/index')
def index():

    return render_template('index.html',
                           title='Home')
