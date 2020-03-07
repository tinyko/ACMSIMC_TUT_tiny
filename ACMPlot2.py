import plotly.graph_objs as go
import plotly.offline as py
from plotly.subplots import make_subplots
import pandas as pd

import numpy as np


def draw_plotly(RUN_TIME, MACHINE_TS_INVERSE):
    df = pd.read_csv("algorithm.dat", sep=",")

    time = np.linspace(0, RUN_TIME, int(RUN_TIME * MACHINE_TS_INVERSE / 2))

    df["time"] = time

    color1 = "#9467bd"
    color2 = "#F08B00"

    # trace1 = go.Scatter(
    #     x=df["time"], y=df["ACM.x0"], name="ACM.x0", line=dict(color=color1)
    # )
    # trace2 = go.Scatter(
    #     x=df["time"], y=df["ACM.x1"], name="ACM.x1", line=dict(color=color2)
    # )
    # data = [trace1, trace2]

    # # 2nd figure
    # trace3 = go.Scatter(
    #     x=df["time"], y=df["ACM.x2"], name="ACM.x2", line=dict(color=color1)
    # )
    # trace4 = go.Scatter(
    #     x=df["time"], y=df["ACM.x3"], name="ACM.x3", line=dict(color=color2)
    # )
    # data2 = [trace3, trace4]

    df["ACM.rpm_mes"] = df["ACM.x4"] * 60.0 / (2 * np.pi * 2)
    trace5 = go.Scatter(
        x=df["time"], y=df["ACM.rpm_cmd"], name="ACM.rpm_cmd", line=dict(color=color1)
    )
    trace6 = go.Scatter(
        x=df["time"], y=df["ACM.rpm_mes"], name="ACM.rpm_mes", line=dict(color=color2)
    )
    data3 = [trace5, trace6]

    trace7 = go.Scatter(
        x=df["time"], y=df["ACM.Tem"], name="ACM.Tem", line=dict(color=color1)
    )

    trace8 = go.Scatter(
        x=df["time"], y=df["ACM.ual"], name="ACM.ual", line=dict(color=color1)
    )
    trace9 = go.Scatter(
        x=df["time"], y=df["ACM.ube"], name="ACM.ube", line=dict(color=color2)
    )
    data4 = [trace8, trace9]

    trace10 = go.Scatter(
        x=df["time"], y=df["ACM.ial"], name="ACM.ial", line=dict(color=color1)
    )
    trace11 = go.Scatter(
        x=df["time"], y=df["ACM.ibe"], name="ACM.ibe", line=dict(color=color2)
    )
    data5 = [trace10, trace11]
    # py.iplot(fig2)
    fig = make_subplots(
        rows=4,
        subplot_titles=(
            r"$Voltage\ in\ \alpha\beta\ frame$",
            r"$Current\ in\ \alpha\beta\ frame$",
            r"$Speed\ Command\ and\ Measurement$",
            r"$Eletromagnetic\ Torque$",
            # r"$error rpm$",
        ),
    )
    trace12 = go.Scatter(
        x=df["time"], y=df["e_omega"], name="e_rpm", line=dict(color=color2)
    )
    fig.add_traces(data=data4, rows=[1, 1], cols=[1, 1])
    fig.add_traces(data=data5, rows=[2, 2], cols=[1, 1])
    fig.add_traces(data=data3, rows=[3, 3], cols=[1, 1])
    fig.add_trace(trace7, row=4, col=1)
    # fig.add_trace(trace12, row=3, col=1)

    fig.update_layout(height=900, showlegend=False)
    fig.show()
