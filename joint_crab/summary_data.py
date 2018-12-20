"""Print a data summary text file."""
import logging
from pathlib import Path
from gammapy.data import DataStore
import astropy.units as u
from .conf import config

log = logging.getLogger(__name__)


def make_summary_data():
    lines = make_summary_data_lines()

    filename = "results/summary/data.md"
    log.info(f"Writing {filename}")
    Path(filename).write_text("\n".join(lines))


def make_summary_data_lines():
    yield "# Data summary\n"
    yield ""

    for name in config.all_datasets:
        for line in make_summary_for_dataset(name):
            yield line


def make_summary_for_dataset(name):
    dataset = config.datasets[name]
    spec = dataset.get_SpectrumObservationList().stack()
    stats = spec.total_stats

    yield f"## {dataset.name}"
    yield "\nInfos from dataset:\n"

    energy_min = dataset.energy_range[0].to("TeV").value
    yield f"- energy_min: {energy_min:.5f} TeV"

    energy_max = dataset.energy_range[1].to("TeV").value
    yield f"- energy_max: {energy_max:.5f} TeV"

    r_on = dataset.on_radius.to("deg").value
    yield f"- r_on: {r_on:.5f} deg"

    yield "\nInfos from spec:\n"
    livetime = spec.livetime.to("hour").value
    yield f"- livetime: {livetime:.2f} hours"

    # n_on = spec.on_vector.data.data.value.sum()
    yield f"- n_on: {stats.n_on}"
    yield f"- excess: {stats.excess:.1f}"
    yield f"- background: {stats.background:.1f}"

    if dataset.name != "fermi":
        data_store = DataStore.from_dir(f"data/{dataset.name}")
        yield "\nInfos from datastore:\n"

        n_obs = len(dataset.obs_ids)
        yield f"- n_obs: {n_obs}"

        obstime = data_store.obs_table["ONTIME"].sum() / 3600
        yield f"- obstime: {obstime:.2f} hours"

    yield ""


def make_summary_latex():
    """Make summary LaTeX table for the paper."""
    path = "results/summary/data.tex"
    log.info(f"Writing {path}")

    table = open(path, "w")
    table.write(r"\begin{tabular}{lrrrrrl}")
    table.write("\n")
    table.write(r"\hline")
    table.write("\n")
    table.write(
        r"Dataset &  $T_{\rm obs}$ & $E_{\rm min}$ & $E_{\rm max}$ & $N_{\rm on}$ & $N_{\rm bkg}$ & $R_{\rm on}$  \\"
    )
    table.write("\n")
    table.write(
        r"        &       & TeV           & TeV           &              &         & deg        \\ \hline"
    )
    table.write("\n")

    for name in config.all_datasets:
        dataset = config.datasets[name]
        e_min_list = []
        for _spec in dataset.get_SpectrumObservationList():
            e_min_list.append(_spec.lo_threshold)
        spec = dataset.get_SpectrumObservationList().stack()
        stats = spec.total_stats

        row_name = rf"\{name} & "
        if name == "fermi":
            T_obs = r"$\sim$7 yr & "
            e_min = dataset.energy_range[0].to("TeV").value
            E_min = f"{e_min:.2f} & "
        else:
            data_store = DataStore.from_dir(f"data/{dataset.name}")
            ontime = sum(data_store.obs_table["ONTIME"]) * u.s
            ontime = ontime.to('h').value
            T_obs = f"{ontime:.2f} h & "
            # in case of the IACT e_min is taken from the staked obs
            e_min = min(e_min_list)
            E_min = f"{e_min.to('TeV').value:.2f} & "

        e_max = dataset.energy_range[1].to("TeV").value
        E_max = f"{e_max:.0f} & "
        N_on = f"{stats.n_on} & "
        N_bkg = f"{stats.background:.1f} & "
        r_on = dataset.on_radius.to("deg").value
        R_on = rf"{r_on:.2f} \\"

        table.write(row_name + T_obs + E_min + E_max + N_on + N_bkg + R_on)
        table.write("\n")

    table.write("\hline")
    table.write("\n")
    table.write("\end{tabular}")
    table.close()
