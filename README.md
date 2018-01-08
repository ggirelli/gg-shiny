gg-shiny
===

Web-interface apps, powered by [R-shiny](https://shiny.rstudio.com/) and [plotly](https://plot.ly/), to easily generate plots for a variety of projects.

### Installation
```
git clone https://github.com/ggirelli/gg-shiny
cd gg-shiny
echo "export PATH:$(pwd)" | tee -a ~/.bash_profile
```

### Usage

First, create an app folder (e.g., `AppName`) with the proper R-shiny structure and files inside. Then, run it with `runApp AppName`. Additionaly, specify port and host Ip with the `--port` and `--host` arguments, separately.