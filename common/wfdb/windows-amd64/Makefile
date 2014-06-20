#Makefile for downloading and building all of the native libraries to support WFDB;
#WFDB will be build on the current directory under the tmp subfolder

#This Makefile will build the following packages:
# curl 
# wfdb (with curl dependency)
# edr  (with wfdb dependency)
# mse
# ecgpuwave
#
# The file will fetch the required files from their respective sources and comp
# and compile everything accordingly. It may be necessary to run this file as 
# with root priveleges.
#
# You will need to have wget and zip installed in Mingw for this Makefile to work
#
#

#Package names and versions
CURL=curl-7.33.0
WFDB="wfdb-10.5.22"
ECGPUWAVE=ecgpuwave-1.3.2

INSTALL_DIR:=$(shell pwd)
BUILD_DIR=build

clean: clean-wfdb clean-curl clean-edr clean-ecgpuwave clean-mse
	rm -rf ./tmp; \
	rm -rf ./$(BUILD_DIR) ;\
	rm -rf ./bin/*	;\
	rm -rf ./lib/*

##CURL TARGET
clean-curl: 
	rm -rf ./tmp/$(CURL); \
	rm -f ./curl.isconfig; \
	rm -f ./curl.isbuilt; \
	rm -rf ./$(BUILD_DIR) ;\
	rm -rf ./tmp; \
	rm -f curl.tar.gz
curl.tar.gz:
	wget http://curl.haxx.se/download/$(CURL).tar.gz; \
	mv $(CURL).tar.gz curl.tar.gz; \
	mkdir -p ./tmp/ ;\
	mkdir ./$(BUILD_DIR) ;\
	cp ./curl.tar.gz ./tmp/;
curl.isconfig: curl.tar.gz
	cd ./tmp/; \
	tar xfvz curl.tar.gz; \
	cd $(CURL); \
	./configure --prefix=$(INSTALL_DIR)/$(BUILD_DIR) --without-ssl --without-zlib --without-libssh2 --disable-ldap; \
	touch ../../curl.isconfig
curl.isbuilt: curl.isconfig
	cd ./tmp/$(CURL); \
	make ;\
	make install ;\
	touch ../../curl.isbuilt;
curl: curl.isbuilt

##WFDB TARGET
#TODO: remove for build directory as well??
clean-wfdb: 
	rm -rf ./tmp/$(WFDB) ;\
	rm -f ./wfdb.isconfig; \
	rm -f ./wfdb.isbuilt; \
	rm -f wfdb.tar.gz
wfdb.tar.gz:
	wget http://www.physionet.org/physiotools/wfdb.tar.gz; \
	cp ./wfdb.tar.gz ./tmp/wfdb.tar.gz;
wfdb.isconfig: wfdb.tar.gz
	cd ./tmp; \
	tar xvfz wfdb.tar.gz; \
	cd $(WFDB); \
	./configure --prefix=$(INSTALL_DIR)/$(BUILD_DIR) --with-libcurl ;\
	touch ../../wfdb.isconfig
wfdb.isbuilt: PATH:=$(PATH):$(INSTALL_DIR)/$(BUILD_DIR)/bin
wfdb.isbuilt: curl.isbuilt wfdb.isconfig
	cd ./tmp/;\
	export LIBRARY_PATH=$(INSTALL_DIR)/$(BUILD_DIR)/lib ;\
	cd ./$(WFDB)/ ;\
	export PATH; \
	export C_INCLUDE_PATH=$C_INCLUDE_PATH:$(INSTALL_DIR)/$(BUILD_DIR)/include; \
	make install ;\
	touch ../../wfdb.isbuilt
wfdb: curl wfdb.isbuilt

##EDR TARGET
clean-edr:
	rm -f edr*; \
	rm -f ./$(BUILD_DIR)/edr
edr: LL=$(INSTALL_DIR)/$(BUILD_DIR)
edr: wfdb edr.c
	mkdir -p $(LL)/include/wfdb ;\
	cp -r $(INSTALL_DIR)/tmp/$(WFDB)/lib/*.h* $(LL)/include/wfdb/; \
	cp -r $(INSTALL_DIR)/tmp/$(WFDB)/lib/*.dll $(LL)/lib/wfdb.dll; \
	gcc -I$(LL)/include/ -o edr -O edr.c -L$(LL)/lib -L$(LL)/lib64 -lm -lwfdb -lcurl
	cp edr ./$(BUILD_DIR)/bin
edr.c:
	wget http://www.physionet.org/physiotools/edr/edr.c


#MULTISCALE ENTROPY TARGET
clean-mse:
	rm -f mse*; \
	rm -f ./$(BUILD_DIR)/mse
mse: curl mse.c
	gcc -o mse -O mse.c -lm ;\
	cp mse ./$(BUILD_DIR)/bin
mse.c:
	wget http://www.physionet.org/physiotools/mse/mse.c
	
	
#ECGPU Wave Target
clean-ecgpuwave: 
	rm -rf ./tmp/$(ECGPUWAVE) ;\
	rm -f ./tmp/ecgpuwave* ;\
	rm -f ./ecgpuwave.isconfig; \
	rm -f ./ecgpuwave.isbuilt; \
	rm -f ecgpuwave*
ecgpuwave.tar.gz:
	wget http://www.physionet.org/physiotools/ecgpuwave/src/$(ECGPUWAVE).tar.gz ;\
	mv $(ECGPUWAVE).tar.gz ecgpuwave.tar.gz ;\
	cp ecgpuwave.tar.gz ./tmp/ecgpuwave.tar.gz ;\
	cd ./tmp; \
	tar xvfz ecgpuwave.tar.gz
ecgpuwave-gfortran-addon.zip: ecgpuwave.tar.gz
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/ecgpuwave-gfortran-addon.zip ;\
	cp ecgpuwave-gfortran-addon.zip ./tmp/ecgpuwave-gfortran-addon.zip ;\
	cd ./tmp; \
	unzip -o ecgpuwave-gfortran-addon.zip ;\
	cp -Lfv ./*.f ./$(ECGPUWAVE)/	
ecgpuwave.isbuilt: C_INCLUDE_PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/include
ecgpuwave.isbuilt: LIBRARY_PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/lib 
ecgpuwave.isbuilt: ecgpuwave-gfortran-addon.zip
	cd ./tmp/$(ECGPUWAVE);\
	export LIBRARY_PATH ;\
	export C_INCLUDE_PATH;\
	export BUILD_DIR ;\
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/ecgpuwave-mingw-Makefile ;\
	mv ecgpuwave-mingw-Makefile Makefile ;\
	make ;\
	cp -v ecgpuwave.exe $(INSTALL_DIR)/$(BUILD_DIR)/bin;\
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/libgfortran-3.dll;\
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/libgcc_s_dw2-1.dll;\
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/libquadmath-0.dll;\
	cp -v libgfortran-3.dll $(INSTALL_DIR)/$(BUILD_DIR)/bin; \
	cp -v libgcc_s_dw2-1.dll $(INSTALL_DIR)/$(BUILD_DIR)/bin ;\
	cp -v libquadmath-0.dll $(INSTALL_DIR)/$(BUILD_DIR)/bin ;\
	touch ../../ecgpuwave.isbuilt
ecgpuwave: wfdb ecgpuwave.isbuilt

all: curl wfdb edr mse ecgpuwave

checkin:
	rm -fv $(INSTALL_DIR)/$(BUILD_DIR)/bin/hr* ;\
	rm -fv $(INSTALL_DIR)/$(BUILD_DIR)/bin/pschart ;\
	rm -fv $(INSTALL_DIR)/$(BUILD_DIR)/bin/psdf ;\
	rm -rfv $(INSTALL_DIR)/$(BUILD_DIR)/lib/pkgconfig/libcurl.pc 

all: curl wfdb mse edr ecgpuwave	

#The install target will copy the contents of the build directory
#to the current lib lib64 and bin directories and remove configuration specific files
install: all
	rm -rv ./bin/* ;\
	rm -rf ./lib/* ;\
	cp -r ./$(BUILD_DIR)/bin/* ./bin/ ;\
	cp -r ./$(BUILD_DIR)/lib/* ./lib/ ;\
	rm -f ./bin/psfd* ./bin/hrmem* ./bin/cur-config* ./bin/hrfft* ;\
	rm -f ./bin/pschar* ./bin/hrlomb* ./bin/hrplot* ./bin/setwfdb* ;\
	rm -f ./bin/cshsetwfdb* ./bin/curl-config* ./lib/libcurl.la* ;\
	rm -f ./lib/pkgconfig/libcurl.pc*

#Package will clean all temporary files 
#generated by this Makefile so that the Toolbox binaries 
#are ready to be deployed
package: install
	rm -rf ./tmp/ ;\
	rm -rf ./$(BUILD_DIR) ;\
	rm -f curl* ecgpuwave* edr* mse* wfdb* *~
	
	