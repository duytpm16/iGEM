# The shinyswatch package provides themes from https://bootswatch.com/

from faicons import icon_svg
import shinyswatch
from pathlib import Path
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import FileInfo
from datatable import dt, f, fread
from plotnine import aes, geom_point, ggplot, scale_x_continuous, scale_y_continuous, scale_color_manual, theme, element_blank, element_line, element_text
import pandas as pd

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
    
    @reactive.Effect
    @reactive.event(input.file)
    def _():
        fi: list[FileInfo] = input.file()
        df = fread(file=fi[0]["datapath"], sep = "\t", header = True)

        df.colindex("P_Value_Joint")
        df[:,f.P_Value_Joint] = df[:,dt.math.log10(f.P_Value_Joint) * -1.0]
        df = df[:,:,dt.sort(f.CHR, f.POS)]
        df['cumulative_pos'] = df[:,dt.cumsum(f.POS)]
        mp = df[:, (dt.max(f.cumulative_pos) + dt.min(f.cumulative_pos)) / 2, dt.by("CHR")]
        mp = mp.to_pandas()
        df = df.to_pandas()
        mid_points.set(mp)
        data.set(df)


    
    @reactive.Effect
    @reactive.event(input.gwis_button)
    def _():
        if not gwis_button_clicked():
            gwis_button_clicked.set(True)
            ui.insert_ui(
                ui.page_fluid(
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
                        ui.row(
                            ui.column(
                                8,
                                ui.div(
                                    {"class": "card mb-4 box-shadow"},
                                    ui.div(
                                        {"class": "card-body", "style": "border-top: 3px solid black;"},
                                        ui.h5({"class": "card-title mt-0"}, "Manhattan Plot"),
                                        ui.output_plot("mh_plot", height="400px", click=True, hover=True, brush=True),
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
                    ),
                ),
                selector="#gwis_button",
                where="afterEnd",
            )



    @output
    @render.plot
    def mh_plot():
        return (
            (
                ggplot(data(), aes("cumulative_pos", "P_Value_Joint")) +
                    geom_point(aes(color='factor(CHR)')) +
                    scale_x_continuous(expand = (0.01,0), 
                                       breaks = mid_points()['C0'], 
                                       labels = list(mid_points()['CHR'])) +
                    scale_color_manual(values = ["black", "grey"]*11, 
                                       guide  = None) +
                    scale_y_continuous(expand = (0.01, 0)) +
                    theme(panel_background=element_blank(),
                          panel_grid=element_line(color="#f7f7f7"),
                          axis_line=element_line(size=0.6),
                          axis_text=element_text(size = 8),
                          axis_title=element_text(size=10))
            )
        )


    @reactive.Effect
    @reactive.event(input.mh_plot_click)
    def _():
        ui.insert_ui(
            ui.page_fluid(
                ui.output_table("test_table")
            ),
            selector="#mh_plot",
            where="afterEnd",
        )

    @output
    @render.table
    def test_table():
        return(mid_points().head(1000))

www_dir = Path(__file__).parent / "www"
app = App(ui=app_ui, server=server, static_assets=www_dir)