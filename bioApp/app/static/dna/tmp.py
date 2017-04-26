import numpy as np
from Bio import SeqIO
import plotly.plotly as py
import plotly.graph_objs as go
from collections import defaultdict



channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12']

seq_file ='A1.fA1_001.ab1'
record = SeqIO.read(seq_file, 'abi')

trace = defaultdict(list)

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
#print(trim_sequence_record.seq)

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

    title = 'A1.fA1_001.ab1 Raw Plots',
    hovermode= 'closest',
    showlegend = True,
    autosize = True,

    font = dict(
        size = 22,
    ),

    xaxis = dict(
        title = 'sequence data',
        range = [0,220],
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
        range = [0,2500],
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

py.plot(fig, filename='ABI Sequencing Test')
