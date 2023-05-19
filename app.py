# The shinyswatch package provides themes from https://bootswatch.com/

from faicons import icon_svg
import shinyswatch
from pathlib import Path
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import FileInfo
from shiny.plotutils import near_points
from datatable import dt, f, fread
from plotnine import aes, geom_point, ggplot, scale_x_continuous, scale_y_continuous, scale_color_manual, theme, element_blank, element_line, element_text
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
from itables import init_notebook_mode
init_notebook_mode(all_interactive=True)

app_ui = ui.page_navbar(
    # Available themes:
    #  cerulean, cosmo, cyborg, darkly, flatly, journal, litera, lumen, lux,
    #  materia, minty, morph, pulse, quartz, sandstone, simplex, sketchy, slate,
    #  solar, spacelab, superhero, united, vapor, yeti, zephyr
    shinyswatch.theme.spacelab(),

    ui.nav(
        "v1.0.0",
        ui.tags.head(
            ui.tags.link(rel="stylesheet", href="radar_style.css"),
        ),
        ui.layout_sidebar(
            ui.panel_sidebar(
                ui.input_file("file", "File input:"),
                width=2,
            ),
            ui.panel_main(
                
                ui.page_fluid(
                    ui.row(
                        ui.column(
                            12, 
                            ui.input_action_button("gwis_button", "GWIS", icon=icon_svg("chart-column",width="30px"), width="25%", class_="btn-info")
                        ),
                    ),
                ),
                ui.HTML("<br/>"),
                # ui.HTML("<div style='height: 200px; overflow: auto;'>"),
                # ui.output_table("table"),
                # ui.HTML("</div>"),
                width=10,
            ),
        ),
    ),
    title=ui.h1("iGEM"),
)


def server(input: Inputs, output: Outputs, session: Session):
    gwis_button_clicked = reactive.Value(False)
    data = reactive.Value()
    mid_points = reactive.Value()
    mh_pt_colors = reactive.Value()
    
    @reactive.Effect
    @reactive.event(input.file)
    def _():
        fi: list[FileInfo] = input.file()
        df = fread(file=fi[0]["datapath"], sep = "\t", header = True)

        df.colindex("P_Value_Joint")
        df[:,f.P_Value_Joint] = df[:,dt.math.log10(f.P_Value_Joint) * -1.0]
        df = df.to_pandas()
        df['POS'] = df['POS'].astype(np.int64)

        running_pos = 0
        cumulative_pos = []
        for chrom, group_df in df.groupby('CHR'):  
            cumulative_pos.append(group_df['POS'] + running_pos)
            running_pos += group_df['POS'].max()
        df['cumulative_pos'] = pd.concat(cumulative_pos)

        df = dt.Frame(df)
        mp = df[:, dt.median(f.cumulative_pos), dt.by("CHR")]
        ct = df[:, dt.count(f.CHR), dt.by("CHR")]
        ct = ct.to_pandas()
        c = []
        color = ["black", "grey"]
        colors = cycle(color)
        for i in ct['CHR.0']:
            color = next(colors)
            c.extend([color for j in range(i)])

        df = df.to_pandas()
        mp = mp.to_pandas()
        mh_pt_colors.set(c)
        
        mid_points.set(mp)
        data.set(df)


    
    @reactive.Effect
    @reactive.event(input.gwis_button)
    def _():
        if not gwis_button_clicked():
            gwis_button_clicked.set(True)
            ui.insert_ui(
                ui.page_fluid(
                    ui.div({"style": "overflow-x: auto;"}),
                    ui.input_select(
                        "gwis_select",
                        "Choose a test:",
                        ["Joint", "Interaction", "Marginal"],
                    ),
                    ui.input_select(
                        "se_select",
                        "Choose standard errors:",
                        ["Model-based", "Robust"],
                    ),
                    ui.panel_conditional(
                        'input.gwis_select === "Joint" && input.se_select === "Model-based"',
                        ui.HTML("<div style='overflow-y: auto; max-height: 575px'>"),
                        ui.row(
                            ui.column(
                                8,
                                ui.div(
                                    {"class": "card mb-4 box-shadow"},
                                    ui.div(
                                        {"class": "card-body", "style": "border-top: 3px solid black;"},
                                        ui.h5({"class": "card-title mt-0"}, "Manhattan Plot"),
                                        ui.output_plot("mh_plot", height="400px", click=True),
                                    ),
                                ),
                            ),
                            ui.column(
                                4,
                                ui.div(
                                    {"class": "card mb-4 box-shadow"},
                                    ui.div(
                                        {"class": "card-body", "style": "border-top: 3px solid black;"},
                                        ui.h5({"class": "card-title mt-0"}, "Quantile-Quantile Plot"),
                                        ui.output_plot("qq_plot", height="400px"),
                                    ),
                                ),
                            ),
                        ),
                        ui.row(
                            ui.column(
                                8,
                                ui.div(
                                    {"class": "card mb-4 box-shadow"},
                                    ui.div(
                                        {"class": "card-body", "style": "border-top: 3px solid black;"},
                                        ui.h5({"class": "card-title mt-0"}, "Variants in Region"),
                                        ui.HTML("<br>"),
                                        ui.HTML("<div style='height: 400px; overflow: scroll;'>"),
                                        ui.output_table("test_table", height="800px"),
                                        ui.HTML("</div>"),
                                    ),
                                ),

                            ),              
                        ),
                    ),
                ),
                selector="#gwis_button",
                where="afterEnd",
            )


    @output
    @render.plot
    def mh_plot():

        _, ax = plt.subplots(figsize=(9, 3))
        ax.scatter(data()["cumulative_pos"], data()["P_Value_Joint"], c=mh_pt_colors(), alpha = 0.6)
        ax.set_xticks(mid_points()["cumulative_pos"])
        ax.set_xticklabels(mid_points()["CHR"])

        ax.set_xlabel("Chromosome")
        # ax.set_ylabel(ylabel)
        ax.margins(x = 0.01)
        return ax

        # return (
        #     (
        #         ggplot(data(), aes("cumulative_pos", "P_Value_Joint")) +
        #             geom_point(aes(color='factor(CHR)')) +
        #             scale_x_continuous(expand = (0.01,0), 
        #                                breaks = mid_points()['cumulative_pos'], 
        #                                labels = list(mid_points()['CHR'])) +
        #             scale_y_continuous(expand = (0.01, 0)) +
        #             scale_color_manual(values = ["black", "grey"]*11, 
        #                                guide  = None) +
        #             theme(panel_background = element_blank(),
        #                   panel_grid = element_line(color = "#f7f7f7"),
        #                   axis_line  = element_line(size  = 0.6),
        #                   axis_text  = element_text(size  = 8),
        #                   axis_title = element_text(size  = 10))
        #     )
        # )


    @output
    @render.table
    def test_table():
        d = near_points(
            data(),
            input.mh_plot_click(),
            xvar="cumulative_pos",
            yvar="P_Value_Joint",
            threshold=5
        )
        
        return near_points(
            data(),
            input.mh_plot_click(),
            xvar="cumulative_pos",
            yvar="P_Value_Joint",
            threshold=5
        )[["SNPID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF"]]

www_dir = Path(__file__).parent / "www"
app = App(ui=app_ui, server=server, static_assets=www_dir)