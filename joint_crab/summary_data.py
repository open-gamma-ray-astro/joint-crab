"""Print a data summary text file."""
import logging
from pathlib import Path
from .conf import config

log = logging.getLogger(__name__)


def make_summary_data():
    lines = make_summary_data_lines()

    filename = "results/summary_data.md"
    log.info(f"Writing {filename}")
    Path(filename).write_text("\n".join(lines))


def make_summary_data_lines():
    yield "# Data summary\n"
    yield ""

    for name in config.all_datasets:
        for line in make_summary_for_dataset(name):
            yield line


def make_summary_for_dataset(name):
    dataset = config.get_dataset(name)
    datastore = dataset.get_DataStore()
    spec = dataset.get_SpectrumObservationList().stack()
    stats = spec.total_stats

    yield f"## {dataset.name}"
    yield "\nInfos from dataset:\n"

    n_obs = len(dataset.obs_ids)
    yield f"- n_obs: {n_obs}"

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

    if datastore:
        yield "\nInfos from datastore:\n"

        # import IPython; IPython.embed(); 1/0
        obstime = datastore.obs_table["ONTIME"].sum() / 3600
        yield f"- obstime: {obstime:.2f} hours"

    yield ""
