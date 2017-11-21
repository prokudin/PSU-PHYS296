from fitlab.resman import RESMAN
from fitlab.mcsamp import MCSAMP
from fitlab.maxlike import ML
from fitlab.speedtest import SPEEDTEST
from tools.tools import load_config


def run(conf_file, task):
    conf = load_config(conf_file)
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

    return conf


if __name__ == "__main__":
    conf = run("fitlab/inputs/upol.py", 2)
