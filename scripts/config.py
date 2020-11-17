import os


if os.geteuid() == 501:
    data_directory = os.path.expanduser("~/GitHub/slim_recipes/data/")
    scripts_directory = os.path.expanduser("~/GitHub/slim_recipes/scripts/")
    analysis_directory = os.path.expanduser("~/GitHub/slim_recipes/analysis/")

else:
    data_directory = os.path.expanduser("/u/project/ngarud/wrshoema/slim_recipes/data/")
    scripts_directory = os.path.expanduser("/u/project/ngarud/wrshoema/slim_recipes/scripts/")
    analysis_directory = os.path.expanduser("/u/project/ngarud/wrshoema/slim_recipes/analysis/")
