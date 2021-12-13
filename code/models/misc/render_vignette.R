# input
args = commandArgs(trailingOnly=TRUE)

# rendering
rmarkdown::render(here::here(args[1]))

# returning file names created to initiate tracking via dvc
print(args[-1])
