from argparse import Namespace
from resman import RESMAN
from mcsamp import MCSAMP
from maxlike import ML
from speedtest import SPEEDTEST
from tools.tools import load_config


def conf_run(conf, task, noexit=True):
    if noexit:
        conf["screen mode"] = "please don't sys.exit"

    conf["resman"] = RESMAN(conf)

    if task == 0:
        SPEEDTEST(conf).run()
    elif task == 1:
        ML(conf).run_minimize()
    elif task == 2:
        ML(conf).run_leastsq()
    elif task == 3:
        MCSAMP(conf).run_nest()
    elif task == 4:
        MCSAMP(conf).run_imc()
    elif task == 5:
        MCSAMP(conf).analysis()
    elif task == 6:
        MCSAMP(conf).simulation()
    elif task == 7:
        MCSAMP(conf).simulation2()
    else:
        raise ValueError("task must be in range(8)")

    return conf


def run(config,
        task=0,
        runid=0,
        file="",
        list=[],
        reaction="sidis",
        noexit=True
        ):

    args = Namespace(
        config=config,
        task=task,
        runid=runid,
        file=file,
        list=list,
        reaction=reaction
    )

    conf = load_config(args.config)
    conf["args"] = args

    return conf_run(conf, task, noexit=noexit)


if __name__ == "__main__":
    conf = run("inputs/upol_hermes_noevolution.py", task=1)
