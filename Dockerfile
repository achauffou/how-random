# The image provided by the Rocker project builds on Ubuntu 20.04 and provides:
# - R version 4.0.3
# - RStudio server (latest)
# - R tidyverse packages from RStudio Package Manager (2021-02-17 archive)
FROM rocker/tidyverse:4.0.3

# Find out the number of cores available for the installation:
ENV NCPUS=${NCPUS:-1}

# Install additional apt packages (listed in .apt_packages):
COPY .apt_packages /home/rstudio/how-random/.apt_packages
RUN apt-get update && apt-get install -y \
  $(cat /home/rstudio/how-random/.apt_packages | tr '\n' ' ')

# Set curl as default methods to download files in R:
RUN echo 'options(download.file.method = "curl")' >> /usr/local/lib/R/etc/Rprofile.site \
  && echo 'options(download.file.extra = "-L")' >> /usr/local/lib/R/etc/Rprofile.site

# Install CMDSTAN 2.26.1 (2021-02-21):
ENV CMDSTAN_REPO=https://github.com/stan-dev/cmdstan/releases/download/v2.26.1/cmdstan-2.26.1.tar.gz
ENV CMDSTAN=/usr/local/cmdstan
RUN mkdir $CMDSTAN \
  && wget -O cmdstan.tar.gz $CMDSTAN_REPO \
  && tar -xzf cmdstan.tar.gz -C $CMDSTAN --strip-components=1 \
  && cd $CMDSTAN \
  && make -j$NCPUS build \
  && cd / \
  && rm cmdstan.tar.gz

# Install cmdstanr 3.0.0 (2020-12-17):
RUN R -e "devtools::install_github('stan-dev/cmdstanr@v0.4.0')" \
  && echo CMDSTAN=$CMDSTAN >> /usr/local/lib/R/etc/Renviron.site

# Install TeX-Live final 2020 archive using TinyTeX:
ENV CTAN_REPO=ftp://tug.org/historic/systems/texlive/2020/tlnet-final
ENV PATH=$PATH:/usr/local/TinyTeX/bin/x86_64-linux
RUN install2.r --error --skipinstalled -r $CRAN -n $NCPUS tinytex \
  && R -e "library(tinytex); tinytex::install_tinytex( \
    dir = '/usr/local/TinyTeX', version = '2020.12', repo = I('"$CTAN_REPO"'))" \
  && echo "export PATH=$PATH" > /etc/environment \
  && chown -R root:staff /usr/local/TinyTeX \
  && chmod -R ugo+rwx /usr/local/TinyTeX

# Install additional TeX packages (listed in .tex_packages):
COPY .tex_packages /home/rstudio/how-random/.tex_packages
RUN tlmgr update --self \
  && tlmgr install $(cat /home/rstudio/how-random/.tex_packages | tr '\n' ' ')

# Install additional R packages (listed in .r_packages):
COPY .r_packages /home/rstudio/how-random/.r_packages
RUN install2.r --error --skipinstalled -r $CRAN -n $NCPUS \
  $(cat /home/rstudio/how-random/.r_packages | tr '\n' ' ')
