# How random ecological interactions are?
*Alain Chauffoureaux, under the supervision of Bernat Bramon Mora*

MSc project on the predictability of ecological interactions carried out at ETH Zürich in 2021.

## Installation
Follow the instructions below to install the project container.

### 1. Requirements
This project uses [docker](https://docs.docker.com/) to build an isolated container for the project that enables to reproduce the results on any operating system without requiring manual software installation and package management. All code should run on any system with docker installed. Please check the following requirements below:
* [Docker](https://docs.docker.com/get-docker/) 20.10 (or any compatible version)
* 250 GB of free storage space on disk
* 8 GB of RAM memory

### 2. Clone or download the project repository
There are to alternatives available to get the code repository:
* **Clone the repository:** [Git](https://git-scm.com/) users can clone the repository with `git clone https://github.com/achauffou/how-random.git`.
* **Download and extract the repository:** Alternatively, download the [latest version](https://github.com/achauffou/how-random) or [any release](https://github.com/achauffou/how-random/releases) of the repository and extract the archive. Do not use this option if you want to contribute to the project.

### 3. Build the docker image
Prior to using the project, it is necessary to build a docker image to contain it. Please be patient, depending on your machine OS and internet connection, it might take sevral minutes to hours to build the image. The image takes up about 4 GB of disk storage. Use the following command from project root folder to build the docker image:
```
docker build -t how-random-image .
```
*Note:* Linux users might need sudo privileges to use docker.

### 4. Run the docker container
Any code from the project is meant to be run within a container based on the previously build docker image. To run a new docker container for the project, use the following command from the project root folder:
```
docker run -d -e DISABLE_AUTH=true -p 8787:8787 -v ${PWD}:/home/rstudio/how-random --name how-random-container how-random-image
```
*Note:* Linux users might need sudo privileges to use docker.

To stop or restart the container, use `docker stop how-random-container` and `docker restart how-random-container`.

### 5. Store NCBI and GBIF creditentials
This project accesses several databases, including GBIF and NCBI that both require authentication information. The user personal NCBI key and GBIF authentication information must be stored in the `.Renviron` file as follows:
```
ENTREZ_KEY='<your_ncbi_api_key>'
GBIF_USER='<your_gbif_username>'
GBIF_PWD='<your_gbif_password>'
GBIF_EMAIL='<your_gbif_email>'
```
To get these creditentials, create free [NCBI](https://www.ncbi.nlm.nih.gov/account/) and [GBIF](https://www.gbif.org/user/profile) user accounts. Once the NCBI user account has been created, follow [these instructions](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) to request a NCBI ENTREZ API key.

## Contribution
### Accessing Rstudio server
Once the project docker container is running, it is possible to develop and run R code directly from Rstudio server. It can be accessed via port 8787 of the local host: [localhost:8787](localhost:8787). To access it on another port, modify the `docker run` command port mapping argument accordingly.

### Opening an interactive terminal with sudo privileges
The terminal from Rstudio server only grants rights for user *rstudio*. If you need to open an interactive terminal in the container with sudo privileges, use the following command:
```
docker exec -it -u root how-random-container bash
```
*Note:* Use `exit` to close the interactive terminal.
