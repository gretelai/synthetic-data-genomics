# This script installs the latest R from the R-project
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
sudo apt update
sudo apt install r-base

# Fix broken dependency in Ubuntu 18.04 by installing libreadline directly
wget http://ftp.debian.org/debian/pool/main/r/readline6/libreadline6_6.3-8+b3_amd64.deb
wget http://ftp.debian.org/debian/pool/main/g/glibc/multiarch-support_2.19-18+deb8u10_amd64.deb
sudo apt install ./libreadline6_6.3-8+b3_amd64.deb ./multiarch-support_2.19-18+deb8u10_amd64.deb 
echo "Cleaning up downloaded packages"
rm -f libreadline*.deb* multiarch*.deb*
