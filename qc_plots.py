import argparse
import os

# import cf_xarray
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# from matplotlib.figure import Figure

# Suppress warning -- writing a lot of figures
plt.rcParams["figure.max_open_warning"] = 0


def plot_results(data, variable_name, results, title, test_name):
    time = data["time"]
    obs = data[variable_name]
    qc_test = results[test_name]

    qc_pass = np.ma.masked_where(qc_test != 1, obs)
    qc_notrun = np.ma.masked_where(qc_test != 2, obs)
    qc_suspect = np.ma.masked_where(qc_test != 3, obs)
    qc_fail = np.ma.masked_where(qc_test != 4, obs)

    fig, ax = plt.subplots(figsize=(15, 3.75))

    ax.set_xlabel("time")
    ax.set_ylabel(f"{data[variable_name].long_name} [{data[variable_name].units}]")

    kw = {"marker": "o", "linestyle": "none"}
    ax.plot(time, obs, label="obs", color="#A6CEE3")
    ax.plot(
        time,
        qc_notrun,
        markersize=2,
        label="qc not run",
        color="gray",
        alpha=0.2,
        **kw,
    )
    ax.plot(
        time,
        qc_pass,
        markersize=4,
        label="qc pass",
        color="green",
        alpha=0.5,
        **kw,
    )
    ax.plot(
        time,
        qc_suspect,
        markersize=4,
        label="qc suspect",
        color="orange",
        alpha=0.7,
        **kw,
    )
    ax.plot(time, qc_fail, markersize=6, label="qc fail", color="red", alpha=1.0, **kw)
    ax.legend(loc="best", bbox_to_anchor=(1.12, 0.75))
    ax.grid(True)
    plt.title(title)
    if not os.path.exists(os.path.join(os.getcwd(), "plots")):
        os.mkdir(os.path.join(os.getcwd(), "plots"))
    canvas = FigureCanvas(fig)
    canvas.print_figure(os.path.join(os.getcwd(), "plots", f"{test_name}.png"))


def clparser() -> argparse.ArgumentParser:
    """Create a parser to handle input arguments and displaying.

    a script specific help message.
    """
    desc_msg = """Create plots of ASV CTD cast data with QARTOD flags."""
    parser = argparse.ArgumentParser(description=desc_msg)
    parser.add_argument("ncfile", help="Path to the input NetCDF file.")
    return parser


def main():
    parser = clparser()
    args = parser.parse_args()
    ncfile = args.ncfile

    data = xr.open_dataset(ncfile)
    variables = {}
    for name in data.cf.standard_names:
        if "status_flag" in name:
            variables.update({name.split()[0]: data.cf.standard_names[name]})

    base_ncfile = os.path.basename(ncfile)
    plot_html = os.path.join("plots", f"{base_ncfile}.html")

    with open(plot_html, "w", encoding="utf-8") as f:
        f.write(
            f"""
            <!DOCTYPE html>
            <html>
                <head>
                    <title>ASV CTD QARTOD Plots</title>
                    <meta charset="utf-8">
                    <meta name="viewport" content="width=device-width, initial-scale=1">
                    <script src="https://kit.fontawesome.com/9d6352f212.js" crossorigin="anonymous"></script>
                    <link rel="icon" type="image/x-icon" href="images/favicon.ico" />
                    <link rel="stylesheet" href="assets/whistler.css" type="text/css" />
                    <!-- Bootstrap CSS (v5) -->
                    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.0/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-gH2yIJqKdNHPEq0n4Mqa/HGKIhSkIHeL5AyhkYV8i59U5AR6csBvApHHNl/vI1Bx" crossorigin="anonymous">
                    <link rel="stylesheet" href="https://cdn.rawgit.com/afeld/bootstrap-toc/v1.0.1/dist/bootstrap-toc.min.css"/>
                </head>
                <body data-bs-spy="scroll" data-bs-target="#toc">
                    <div class="container-fluid m-2 p-2">
                        <h1>Plots for {base_ncfile}</h1>
                        <hr>
                        <div class="row">
                        <div class="col-sm-3 p-3">
                            <nav id="toc" data-toggle="toc" class="sticky-top"></nav>
                        </div>
                        <div class="col-sm-9">

        """,
        )
        print("Generating plots:")
        for var in variables:
            print(f"  {var}")
            f.write(f"""<h3>{var}</h3>""")
            for qc in variables[var]:
                f.write(f"""<div class="row"><h5>{qc}</h5>""")
                plot_results(
                    data,
                    var,
                    data,
                    qc,
                    qc,
                )
                f.write(f"""<img src='{qc}.png'>""")
                f.write("""</div>""")
            f.write("<hr>")

        f.write(
            """
                </div>
                </div>
                <!-- jQuery first, then Popper.js + Bootstrap JS -->
                <script src="https://code.jquery.com/jquery-3.5.1.slim.min.js" integrity="sha384-DfXdz2htPH0lsSSs5nCTpuj/zy4C+OGpamoFVy38MVBnE+IbbVYUew+OrCXaRkfj" crossorigin="anonymous"></script>
                <!-- JavaScript Bundle with Popper -->
                <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.2.0/dist/js/bootstrap.bundle.min.js" integrity="sha384-A3rJD856KowSb7dwlZdYEkO39Gagi7vIsF0jrRAoQmDKKtQBHUuLZ9AsSv4jD4Xa" crossorigin="anonymous"></script>
                <script src="https://cdn.rawgit.com/afeld/bootstrap-toc/v1.0.1/dist/bootstrap-toc.min.js"></script>
                <!-- ajax -->
                <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.1.0/jquery.min.js"></script>
                </body>
            </html>
                """,
        )
    print(f"Output created: {plot_html}")


if __name__ == "__main__":
    main()
