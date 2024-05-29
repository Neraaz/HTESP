#!/usr/bin/env python
"""
This script generates an interactive band structure and density of states (DOS) plot
for a given material using Plotly.

Usage: python plot_band_dos_QE.py <material_id> <compound_name> <nkpoint>
"""

import sys
import numpy as np
import plotly.graph_objs as go
from ase.io import espresso
from kpath import kpath

def plot(file, comp):
    """
    Plots the band structure and DOS for the given compound.

    Parameters:
    file (str): Path to the input file for the kpath calculation.
    comp (str): Compound name.
    """
    n = int(sys.argv[3])
    kcut = 0
    pt1, sympt, symlb, pt, _, _ = kpath(file, n, kcut)
    pt = np.array(pt)
    sympt = np.array(sympt)
    
    # Read Fermi energy from the DOS file
    with open(f"{comp}.dos", "r") as f:
        fermi = float(f.readlines()[0].split()[8])
    
    # Load DOS data
    datados = np.loadtxt(f"{comp}.dos", skiprows=1)
    datadosx = datados[:, 1]
    datadosy = datados[:, 0] - fermi

    # Load band structure data
    data = np.loadtxt(f'{comp}.dat.gnu')
    nband = int(data.shape[0] / n)
    print(f"nkpoints: {n} - nband: {nband}")

    # Create band structure traces
    band_traces = []
    for i in range(nband):
        trace = go.Scatter(
            x=pt,
            y=data[0 + n * i:n * (i + 1)][:, 1] - fermi,
            mode='lines',
            line=dict(color='blue'),
            name=f'Band {i+1}'
        )
        band_traces.append(trace)
    
    # Create DOS trace
    dos_trace = go.Scatter(
        x=datadosx,
        y=datadosy,
        mode='lines',
        line=dict(color='orange'),
        name='Total DOS'
    )

    # Create vertical lines at high-symmetry points
    vlines = []
    for sym in sympt:
        vline = go.Scatter(
            x=[sym, sym],
            y=[-5, 5],
            mode='lines',
            line=dict(color='black', dash='dash'),
            showlegend=False
        )
        vlines.append(vline)

    # Layout for the band structure plot
    band_layout = go.Layout(
        title=f'Band Structure for {comp}',
        xaxis=dict(
            title='k-path',
            tickvals=sympt,
            ticktext=symlb
        ),
        yaxis=dict(title='E - Ef (eV)', range=[-5, 5]),
        showlegend=False
    )

    # Layout for the DOS plot
    dos_layout = go.Layout(
        title=f'DOS for {comp}',
        xaxis=dict(title='DOS'),
        yaxis=dict(title='E - Ef (eV)', range=[-5, 5]),
        showlegend=True
    )

    # Combine all traces and create subplots
    fig = go.Figure()
    for trace in band_traces + vlines:
        fig.add_trace(trace)
    fig.add_trace(dos_trace)

    fig.update_layout(
        xaxis_title="k-path / DOS",
        yaxis_title="E - Ef (eV)",
        showlegend=True,
        title=f'Band Structure and DOS for {comp}',
        xaxis=dict(domain=[0, 0.7]),  # Band structure on the left
        xaxis2=dict(domain=[0.75, 1]),  # DOS on the right
        yaxis2=dict(anchor='x2')  # Align y-axis of DOS with its x-axis
    )

    # Update traces to use the appropriate axes
    for i in range(len(band_traces + vlines)):
        fig['data'][i].update(xaxis='x1', yaxis='y1')
    fig['data'][-1].update(xaxis='x2', yaxis='y2')

    # Save plot as HTML file
    fig.write_html(comp + "-band_dos.html")

if __name__ == "__main__":
    print("Usage: python plot_band_dos_QE.py <material_id> <compound_name> <nkpoint>\n")
    mpid = sys.argv[1]
    comp = sys.argv[2]
    filename = 'scf.in'
    plot(filename, comp)
