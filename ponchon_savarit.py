# Install required libraries in Google Colab (if not already installed)
!pip install numpy plotly scipy pandas

import numpy as np
import plotly.graph_objects as go
from scipy.interpolate import interp1d, CubicSpline
from scipy.optimize import root_scalar
import pandas as pd

# Default VLE and enthalpy data (same as in the original HTML code)
default_data = {
    'xData': [0, 0.08, 0.18, 0.25, 0.49, 0.65, 0.79, 0.91, 1.0],
    'yData': [0, 0.28, 0.43, 0.51, 0.73, 0.83, 0.90, 0.96, 1.0],
    'Hl': [24.3, 24.1, 23.2, 22.8, 22.05, 21.75, 21.7, 21.6, 21.4],
    'Hv': [61.2, 59.6, 58.5, 58.1, 56.5, 55.2, 54.4, 53.8, 53.3]
}

# Input parameters (default values from the HTML code)
zF = 0.42  # Feed composition (mole fraction)
F = 10.0   # Feed flow rate (kmol/hr)
xD = 0.97  # Distillate composition (mole fraction)
xW = 0.01  # Bottoms composition (mole fraction)
q = 1.0    # Feed thermal condition (0 to 1)
R = 2.5    # Reflux ratio (dimensionless)

# Interpolation functions
def linear_interpolate(x, x_values, y_values):
    if not x_values or not y_values or len(x_values) != len(y_values):
        return np.nan
    if x <= x_values[0]:
        return y_values[0]
    if x >= x_values[-1]:
        return y_values[-1]
    f = interp1d(x_values, y_values, kind='linear', fill_value="extrapolate")
    return f(x)

def cubic_interpolate(x, x_values, y_values):
    if not x_values or not y_values or len(x_values) != len(y_values) or len(x_values) < 4:
        return np.nan
    if x < x_values[0] or x > x_values[-1]:
        return np.nan
    f = CubicSpline(x_values, y_values, extrapolate=False)
    return f(x)

# Validate input data
def validate_data(data):
    if not all(key in data for key in ['xData', 'yData', 'Hl', 'Hv']):
        return False, "Data must contain xData, yData, Hl, and Hv arrays."
    length = len(data['xData'])
    if (not all(isinstance(arr, (list, np.ndarray)) for arr in [data['xData'], data['yData'], data['Hl'], data['Hv']]) or
        len(data['yData']) != length or len(data['Hl']) != length or len(data['Hv']) != length or length < 4):
        return False, "All arrays must be of equal length and contain at least 4 points."
    for i in range(length):
        if not all(np.isfinite(val) for val in [data['xData'][i], data['yData'][i], data['Hl'][i], data['Hv'][i]]):
            return False, "All values must be finite numbers."
        if data['xData'][i] < 0 or data['xData'][i] > 1 or data['yData'][i] < 0 or data['yData'][i] > 1:
            return False, "xData and yData must be between 0 and 1."
        if i > 0:
            if data['xData'][i] <= data['xData'][i-1]:
                return False, "xData must be monotonically increasing."
            if data['yData'][i] <= data['yData'][i-1]:
                return False, "yData must be monotonically increasing."
    return True, ""

# Stage calculation function
def calculate_stages(xD, xW, zF, HF, q, R, D, W, xDeltaR, HDeltaR, xDeltaS, HDeltaS, data):
    y = xD
    stages = 0
    stage_points = [{'x': xD, 'y': cubic_interpolate(xD, data['yData'], data['Hv'])}]
    tie_lines = []
    construction_lines = []
    stage_compositions = []
    in_rectifying = True
    feed_stage = 0
    error = ""

    while stages < 100:
        # Find x_n from y_n using equilibrium data
        def find_x(x):
            return cubic_interpolate(x, data['xData'], data['yData']) - y
        sol = root_scalar(find_x, bracket=[0, 1], method='bisect')
        if not sol.converged:
            error = f"Stage {stages + 1}: Invalid liquid composition."
            break
        x_n = sol.root

        if x_n <= xW:
            HLxW = linear_interpolate(xW, data['xData'], data['Hl'])
            if np.isfinite(HLxW):
                stage_points.append({'x': xW, 'y': HLxW})
                tie_lines.append({'x': [xW, y], 'y': [HLxW, cubic_interpolate(y, data['yData'], data['Hv'])]})
                stage_compositions.append({'x': xW, 'y': y})
                stages += 1
            break

        HLx_n = linear_interpolate(x_n, data['xData'], data['Hl'])
        if not np.isfinite(HLx_n):
            error = f"Stage {stages + 1}: Failed to interpolate liquid enthalpy."
            break
        stage_points.append({'x': x_n, 'y': HLx_n})

        HVy_n = cubic_interpolate(y, data['yData'], data['Hv'])
        if not np.isfinite(HVy_n):
            error = f"Stage {stages + 1}: Failed to interpolate vapor enthalpy."
            break
        tie_lines.append({'x': [x_n, y], 'y': [HLx_n, HVy_n]})
        stage_compositions.append({'x': x_n, 'y': y})
        stages += 1

        if in_rectifying and x_n <= zF:
            in_rectifying = False
            feed_stage = stages

        xDelta = xDeltaR if in_rectifying else xDeltaS
        HDelta = HDeltaR if in_rectifying else HDeltaS

        if abs(x_n - xDelta) < 1e-6:
            error = f"Stage {stages + 1}: Composition too close to difference point."
            break
        slope = (HDelta - HLx_n) / (xDelta - x_n)
        if not np.isfinite(slope):
            error = f"Stage {stages + 1}: Invalid slope calculation."
            break

        def find_y(y):
            return cubic_interpolate(y, data['yData'], data['Hv']) - (HDelta + slope * (y - xDelta))
        sol = root_scalar(find_y, bracket=[xW, xD], method='bisect')
        if not sol.converged:
            yMin = max(0, y - 0.2)
            yMax = min(1, y + 0.2)
            sol = root_scalar(find_y, bracket=[yMin, yMax], method='bisect')
            if not sol.converged:
                error = f"Stage {stages + 1}: Failed to find valid y_{stages + 2}."
                break
        yNext = sol.root
        HVyNext = cubic_interpolate(yNext, data['yData'], data['Hv'])
        if not np.isfinite(HVyNext):
            error = f"Stage {stages + 1}: Failed to interpolate vapor enthalpy."
            break
        stage_points.append({'x': yNext, 'y': HVyNext})

        if in_rectifying:
            def find_x_end(x):
                return linear_interpolate(x, data['xData'], data['Hl']) - (HDelta + slope * (x - xDelta))
            sol = root_scalar(find_x_end, bracket=[0, 1], method='bisect')
            if not sol.converged:
                error = f"Stage {stages + 1}: Failed to find liquid line intersection."
                break
            xEnd = sol.root
            HEnd = linear_interpolate(xEnd, data['xData'], data['Hl'])
        else:
            xEnd = yNext
            HEnd = HVyNext
        if not np.isfinite(HEnd):
            error = f"Stage {stages + 1}: Invalid construction line endpoint."
            break
        construction_lines.append({'x': [xDelta, xEnd], 'y': [HDelta, HEnd]})
        y = yNext

    if stages == 0 or len(stage_points) < 2:
        error = "Failed to calculate stages. Check input data."
    return stages, stage_points, feed_stage, tie_lines, construction_lines, stage_compositions, error

# Main calculation function
def calculate():
    data = default_data
    valid, error = validate_data(data)
    if not valid:
        print(f"Error: {error}")
        return None, None

    # Validate inputs
    if not all(np.isfinite([zF, F, xD, xW, q, R])):
        print("Error: All inputs must be valid numbers.")
        return None, None
    if F <= 0 or zF < 0 or zF > 1 or xD < 0 or xD > 1 or xW < 0 or xW > 1 or q < 0 or q > 1 or R < 0:
        print("Error: Invalid input values (0 ≤ compositions, q ≤ 1, F > 0, R ≥ 0).")
        return None, None
    if abs(xD - xW) < 1e-6:
        print("Error: Distillate and bottoms compositions must be different.")
        return None, None
    if xD <= zF or zF <= xW:
        print("Error: Compositions must satisfy x_W < z_F < x_D.")
        return None, None

    # Calculate enthalpies
    yF = cubic_interpolate(zF, data['xData'], data['yData'])
    HLzF = linear_interpolate(zF, data['xData'], data['Hl'])
    HVzF = cubic_interpolate(yF, data['yData'], data['Hv'])
    HF = q * HLzF + (1 - q) * HVzF
    HD = linear_interpolate(xD, data['xData'], data['Hl'])
    HW = linear_interpolate(xW, data['xData'], data['Hl'])

    if not all(np.isfinite([yF, HLzF, HVzF, HF, HD, HW])):
        print("Error: Interpolation failed. Ensure compositions are within valid range.")
        return None, None

    # Material balance
    D = F * (zF - xW) / (xD - xW)
    W = F - D
    if not all(np.isfinite([D, W])) or D < 0 or W < 0:
        print("Error: Invalid material balance.")
        return None, None

    # Condenser duty
    HVxD = cubic_interpolate(xD, data['yData'], data['Hv'])
    if not np.isfinite(HVxD):
        print("Error: Failed to calculate vapor enthalpy for x_D.")
        return None, None
    Qc = D * (HVxD - HD) * (R + 1)  # MJ/hr
    QcKW = Qc * 0.27778  # Convert to kW
    xDeltaR = xD
    HDeltaR = HD + Qc / D
    if not all(np.isfinite([Qc, HDeltaR])):
        print("Error: Failed to calculate condenser duty or Δ_R.")
        return None, None

    # Stripping difference point
    xDeltaS = xW
    slope = (HDeltaR - HF) / (xDeltaR - zF)
    HDeltaS = HF + slope * (xDeltaS - zF)
    if not np.isfinite(HDeltaS):
        print("Error: Failed to calculate Δ_S point.")
        return None, None

    # Reboiler duty
    Qr = W * (HW - HDeltaS)  # MJ/hr
    QrKW = Qr * 0.27778  # Convert to kW
    if not np.isfinite(Qr):
        print("Error: Failed to calculate reboiler duty.")
        return None, None

    # Minimum reflux calculations
    yFMin = cubic_interpolate(zF, data['xData'], data['yData'])
    HVyF = cubic_interpolate(yFMin, data['yData'], data['Hv'])
    if not all(np.isfinite([yFMin, HVyF])):
        print("Error: Failed to calculate y_F or H_V(y_F) for minimum reflux.")
        return None, None
    slopeMin = (HVyF - HF) / (yFMin - zF)
    if not np.isfinite(slopeMin):
        print("Error: Failed to calculate slope for minimum reflux line.")
        return None, None
    QPrimeMin = HF + slopeMin * (xD - zF)
    QDoublePrimeMin = HF + slopeMin * (xW - zF)
    if not all(np.isfinite([QPrimeMin, QDoublePrimeMin])):
        print("Error: Failed to calculate minimum reflux intersection points.")
        return None, None
    RMin = (QPrimeMin - HVxD) / (HVxD - HD)
    if not np.isfinite(RMin) or RMin < 0:
        print("Error: Invalid minimum reflux ratio.")
        return None, None

    # Calculate stages
    stages, stage_points, feed_stage, tie_lines, construction_lines, stage_compositions, error = calculate_stages(
        xD, xW, zF, HF, q, R, D, W, xDeltaR, HDeltaR, xDeltaS, HDeltaS, data
    )
    if error:
        print(f"Error: {error}")
        return None, None

    # Plotting
    y_values = data['Hl'] + data['Hv'] + [HF, HD, HW, HDeltaR, HDeltaS, QPrimeMin, QDoublePrimeMin, HVyF]
    y_values = [y for y in y_values if np.isfinite(y)]
    x_values = [zF, xD, xW, xDeltaR, xDeltaS, yFMin] + data['xData'] + data['yData']
    x_values = [x for x in x_values if np.isfinite(x)]
    yMin, yMax = min(y_values) - 10, max(y_values) + 10
    xMin, xMax = min(x_values) - 0.01, max(x_values) + 0.01

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=data['xData'], y=data['Hl'], mode='lines+markers', name='Saturated Liquid',
                             line=dict(color='#1e40af', width=3), marker=dict(size=6),
                             hovertemplate='x: %%{x:.2f}<br>H_L: %%{y:.2f} MJ/kmol'))
    fig.add_trace(go.Scatter(x=data['yData'], y=data['Hv'], mode='lines+markers', name='Saturated Vapor',
                             line=dict(color='#b91c1c', width=3), marker=dict(size=6),
                             hovertemplate='y: %%{x:.2f}<br>H_V: %%{y:.2f} MJ/kmol'))
    fig.add_trace(go.Scatter(x=[xW, xW], y=[yMin, yMax], mode='lines', name='x_W Line',
                             line=dict(color='#475569', width=1.5, dash='dot')))
    fig.add_trace(go.Scatter(x=[xD, xD], y=[yMin, yMax], mode='lines', name='x_D Line',
                             line=dict(color='#475569', width=1.5, dash='dot')))
    fig.add_trace(go.Scatter(x=[xDeltaR, xDeltaS], y=[HDeltaR, HDeltaS], mode='markers+text', name='Difference Points',
                             marker=dict(color='#f59e0b', size=10, symbol='star'),
                             text=['Δ_R', 'Δ_S'], textposition=['middle left', 'middle right'],
                             textfont=dict(size=14, color='#f59e0b'),
                             hovertemplate='%%{x:.2f}, %%{y:.2f} MJ/kmol'))
    fig.add_trace(go.Scatter(x=[zF, xDeltaR, xDeltaS], y=[HF, HDeltaR, HDeltaS], mode='lines', name='Operating Line',
                             line=dict(color='#0891b2', width=2.5),
                             hovertemplate='%%{x:.2f}, %%{y:.2f} MJ/kmol'))
    fig.add_trace(go.Scatter(x=[xW, zF, yFMin, xD], y=[QDoublePrimeMin, HF, HVyF, QPrimeMin], mode='lines', name='Minimum Reflux Line',
                             line=dict(color='#22c55e', width=2, dash='dash'),
                             hovertemplate='%%{x:.2f}, %%{y:.2f} MJ/kmol'))
    fig.add_trace(go.Scatter(x=[xD, xW], y=[QPrimeMin, QDoublePrimeMin], mode='markers+text', name='Min Difference Points',
                             marker=dict(color='#22c55e', size=8, symbol='star-diamond'),
                             text=['Δ_R min', 'Δ_S min'], textposition=['middle left', 'middle right'],
                             textfont=dict(size=14, color='red'),
                             hovertemplate='%%{x:.2f}, %%{y:.2f} MJ/kmol'))

    for i, tie_line in enumerate(tie_lines):
        fig.add_trace(go.Scatter(x=tie_line['x'], y=tie_line['y'], mode='lines', name=f'Stage {i+1}',
                                 line=dict(color='#000000', width=2),
                                 hovertemplate=f'Stage {i+1}<br>(%%{{x:.2f}}, %%{{y:.2f}} MJ/kmol)'))

    for i, const_line in enumerate(construction_lines):
        fig.add_trace(go.Scatter(x=const_line['x'], y=const_line['y'], mode='lines', name=f'Construction Line {i+1}',
                                 line=dict(color='#9333ea', width=1.5, dash='dot'),
                                 hovertemplate=f'Construction Line {i+1}<br>(%%{{x:.2f}}, %%{{y:.2f}} MJ/kmol)'))

    fig.update_layout(
        title='Ponchon-Savarit Diagram for Binary Distillation by ChemEnggCalc',  
        xaxis=dict(title='Mole Fraction (x or y)', range=[xMin, xMax], tickfont=dict(size=12), gridcolor='#e5e7eb'),
        yaxis=dict(title='Enthalpy (MJ/kmol)', range=[yMin, yMax], tickfont=dict(size=12), gridcolor='#e5e7eb'),
        showlegend=False,
        hovermode='closest',
        paper_bgcolor='#ffffff',
        plot_bgcolor='white',
        margin=dict(t=50, b=80, l=80, r=40),
        width=1000,  
        height=600, 
        annotations=[
            dict(
                text='© 2025 chemenggcalc.com | For educational purposes only, not for commercial use.',
                xref='paper', yref='paper', x=0.88, y=0.12, showarrow=False,
                font=dict(size=14, color='rgba(0, 0, 0, 0.5)'),
                xanchor='right', yanchor='bottom'
            )
        ],
        shapes=[dict(type='rect', xref='paper', yref='paper', x0=0, y0=0, x1=1, y1=1,
                     line=dict(color='#1e293b', width=2), fillcolor='rgba(0, 0, 0, 0)')]
    )

    # Display the plot first
    fig.show()

    # Prepare results for table
    results = {
        'Parameter': [
            'Distillate Flow Rate (D)',
            'Bottoms Flow Rate (W)',
            'Rectifying Difference Point (Δ_R)',
            'Stripping Difference Point (Δ_S)',
            'Condenser Duty (Q_C)',
            'Reboiler Duty (Q_R)',
            'Δ_R min',
            'Δ_S min',
            'Minimum Reflux Ratio (R_min)',
            'Number of Theoretical Stages',
            'Feed Stage'
        ],
        'Value': [
            f'{D:.2f} kmol/hr',
            f'{W:.2f} kmol/hr',
            f'({xDeltaR:.2f}, {HDeltaR:.2f} MJ/kmol)',
            f'({xDeltaS:.2f}, {HDeltaS:.2f} MJ/kmol)',
            f'{QcKW:.2f} kW',
            f'{QrKW:.2f} kW',
            f'({xD:.2f}, {QPrimeMin:.2f} MJ/kmol)',
            f'({xW:.2f}, {QDoublePrimeMin:.2f} MJ/kmol)',
            f'{RMin:.2f}',
            f'{stages}',
            f'{feed_stage}'
        ]
    }

    # Create DataFrame for main results
    results_df = pd.DataFrame(results)

    # Create DataFrame for stage compositions
    stage_data = {
        'Stage': [f'Stage {i+1}' for i in range(len(stage_compositions))],
        'x': [f'{stage["x"]:.3f}' for stage in stage_compositions],
        'y': [f'{stage["y"]:.3f}' for stage in stage_compositions]
    }
    stages_df = pd.DataFrame(stage_data)

    # Display tables
    print("\n=== ChemEnggCalc Ponchon-Savarit Calculator Results ===")
    print("\nMain Results:")
    display(results_df)
    print("\nStage Compositions:")
    display(stages_df)

    return fig, {'D': D, 'W': W, 'xDeltaR': xDeltaR, 'HDeltaR': HDeltaR, 'xDeltaS': xDeltaS, 'HDeltaS': HDeltaS,
                 'QcKW': QcKW, 'QrKW': QrKW, 'QPrimeMin': QPrimeMin, 'QDoublePrimeMin': QDoublePrimeMin,
                 'RMin': RMin, 'stages': stages, 'feed_stage': feed_stage, 'stage_compositions': stage_compositions}

# Run the calculation
fig, results = calculate()
