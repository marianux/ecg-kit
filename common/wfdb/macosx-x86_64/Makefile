#Makefile for downloading and building all of the native libraries to support WFDB;
#WFDB will be build on the current directory under the tmp subfolder

#This Makefile will build the following packages:
# curl 
# wfdb (with curl dependency)
# edr  (with wfdb dependency)
# wfdb 32 bit (no curl depency)
# ecgpuwave (with wfdb 32 bit dependency)
# mse
# 
#
# The file will fetch the required files from their resepective sources 
# and compile everything accordingly. 

#This Makefile is not fully functional ! It is provided mainly as a aid for the compilation.
#See comments on respective targets for more information.

#To check linking run: otool -L ./nativelibs/macosx-x86_64/lib/libwfdb.10.dylib

#Package names and versions
CURL=curl-7.33.0
WFDB=wfdb-10.5.22
ECGPUWAVE=ecgpuwave-1.3.2
BUILD_DIR=build
MATLAB=/Applications/MATLAB_R2013b.app/bin/matlab

INSTALL_DIR:=$(shell pwd)
BUILD_DIR32=$(BUILD_DIR)32
WFDB32=$(WFDB)32

clean: clean-wfdb clean-curl clean-edr clean-ecgpuwave
	rm -rf ./tmp; \
	rm -rf ./$(BUILD_DIR);\
	rm -rf ./$(BUILD_DIR) 

##CURL TARGET
clean-curl: 
	rm -rf ./tmp/$(CURL); \
	rm -f ./curl.isconfig; \
	rm -f ./curl.isbuilt; \
	rm -f curl.tar.gz
curl.tar.gz:
	wget http://curl.haxx.se/download/$(CURL).tar.gz; \
	mv $(CURL).tar.gz curl.tar.gz; \
	mkdir -p ./tmp/ ;\
	mkdir -p ./$(BUILD_DIR) ;\
	cp ./curl.tar.gz ./tmp/
curl.isconfig: curl.tar.gz
	cd ./tmp/; \
	tar xfvz curl.tar.gz; \
	cd $(CURL); \
	./configure --prefix=$(INSTALL_DIR)/$(BUILD_DIR) --without-ssl --without-zlib --without-libssh2 --disable-ldap; \
	touch ../../curl.isconfig
curl.isbuilt: curl.isconfig
	cd ./tmp/$(CURL); \
	make ;\
	make install \
	&& touch ../../curl.isbuilt;
curl: curl.isbuilt

##WFDB TARGET
# Requires mannuall editing of the ./conf/darwin*.def files 
# to remove any references to ppc and i386
clean-wfdb: 
	rm -rf ./tmp/$(WFDB)
	rm -f ./wfdb.isconfig; \
	rm -f ./wfdb.isbuilt; \
	rm -f wfdb.tar.gz
wfdb.tar.gz:
	wget http://www.physionet.org/physiotools/wfdb.tar.gz; \
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/darwin-slib.def ;\
        wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/darwin.def ;\
        cp ./wfdb.tar.gz ./tmp/wfdb.tar.gz
wfdb.isconfig: wfdb.tar.gz
	cd ./tmp; \
	tar xvfz wfdb.tar.gz; \
	cp -vf ../darwin.def ./$(WFDB)/conf/ ;\
	cp -vf ../darwin-slib.def ./$(WFDB)/conf/ ;\
	cd $(WFDB); \
	./configure --prefix=$(INSTALL_DIR)/$(BUILD_DIR) --with-libcurl \
	&& touch ../../wfdb.isconfig
wfdb.isbuilt: C_INCLUDE_PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/include
wfdb.isbuilt: LIBRARY_PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/lib
wfdb.isbuilt: PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/bin:$(PATH)
wfdb.isbuilt: curl.isbuilt wfdb.isconfig
	cd ./tmp/;\
	cd ./$(WFDB)/ ;\
	make ;\
	make install \
	&& touch ../../wfdb.isbuilt
wfdb: wfdb.isbuilt

##EDR TARGET
clean-edr:
	rm -f edr*; \
	rm -f ./$(BUILD_DIR)/edr;
edr: LC:=$(shell $(INSTALL_DIR)/$(BUILD_DIR)/bin/curl-config --cflags)
edr: LL=$(INSTALL_DIR)/$(BUILD_DIR)
edr: edr.c wfdb.isbuilt
	gcc $(LC)/curl -o edr -O edr.c -L$(LL)/lib -L$(LL)/lib64 -lm -lwfdb -lcurl
	cp edr ./$(BUILD_DIR)/bin;
edr.c:
	wget http://www.physionet.org/physiotools/edr/edr.c

#MULTISCALE ENTROPY TARGET
clean-mse:
	rm -f mse*; \
	rm -f ./$(BUILD_DIR)/mse;
mse: mse.c
	gcc -o mse -O mse.c -lm
	cp mse ./$(BUILD_DIR)/bin;
mse.c:
	wget http://www.physionet.org/physiotools/mse/mse.c


#WFDB 32 - required by ECGPUWAVE
#To compile i386 version on a Mac, you need to modify the ./conf/darwin-slib.def file
# and pass the -m32 argument to the GCC flag and change all arch types from
# x86_64 to i386 (configure script does not seem to work). After
#the library is build, copy the library files to the buil32 directory. The target 
#below does not work, but is kept here for documentation. TODO: fix this target!
clean-wfdb32:
	rm -rf ./tmp/$(WFDB32)
	rm -f ./wfdb32.isconfig; \
	rm -f ./wfdb32.isbuilt; \
	rm -rf ./$(BUILD_DIR32) ;\
	rm -f wfdb32.tar.gz
wfdb32.isconfig: wfdb.tar.gz
	mkdir -p ./$(BUILD_DIR32) ;\
	cp ./wfdb.tar.gz ./$(BUILD_DIR32)/wfdb.tar.gz ;\
	cd ./$(BUILD_DIR32);\
	tar xvfz wfdb.tar.gz;\
	mv $(WFDB) ../tmp/$(WFDB32) ;\
	cp -vf ../darwin.def ../tmp/$(WFDB32)/conf/ ;\ 
	cp -vf ../darwin-slib.def ../tmp/$(WFDB32)/conf/ ;\
	cd ../tmp/$(WFDB32) ;\
	./configure --prefix=$(INSTALL_DIR)/$(BUILD_DIR32) -m32 --without-netfiles \
	&& touch ../../wfdb32.isconfig
wfdb32.isbuilt: C_INCLUDE_PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/include
wfdb32.isbuilt: LIBRARY_PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/lib
wfdb32.isbuilt: PATH:=$(INSTALL_DIR)/$(BUILD_DIR)/bin:$(PATH)
wfdb32.isbuilt: curl.isbuilt wfdb32.isconfig
	cd ./tmp/$(WFDB32);\
	make install \
	&& touch ../../wfdb32.isbuilt
wfdb32: wfdb32.isbuilt


#ECGPU Wave Target 
#Copy the contents for the BUILD_32/lib to BUILD/lib directory
#For mac compilations you may get gfortran errors like:
# "error trying to exec 'f951': execcvp" . The following link helped resolved this issue:
# http://eftrunk.blogspot.com/2011/11/how-to-fix-error-trying-to-exec-f951.html
#Essentially requiring a re-install of GCC and gfortran tools. 

clean-ecgpuwave: 
	rm -f ecgpuwave.isconfig; \
	rm -f ecgpuwave.isbuilt; \
	rm -f ecgpuwave.tar.gz ;\
	rm -f ecgpuwaveMakefile 
ecgpuwave.tar.gz: 
	wget http://www.physionet.org/physiotools/ecgpuwave/src/$(ECGPUWAVE).tar.gz ;\
	mv $(ECGPUWAVE).tar.gz ecgpuwave.tar.gz ;\
	cp ecgpuwave.tar.gz ./tmp/ecgpuwave.tar.gz;\
	cd ./tmp;\
	tar xvfz ecgpuwave.tar.gz
ecgpuwave-gfortran-addon.zip: ecgpuwave.tar.gz
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/ecgpuwave-gfortran-addon.zip ;\
	wget http://wfdb-app-toolbox.googlecode.com/svn/trunk/dependencies/ecgpuwaveMakefile ;\
	cp ecgpuwave-gfortran-addon.zip ./tmp/ecgpuwave-gfortran-addon.zip ;\
	cp ecgpuwaveMakefile ./tmp/ecgpuwaveMakefile ;\
	cd ./tmp; \
	unzip ecgpuwave-gfortran-addon.zip ;\
	cp -rv ./*.f ./$(ECGPUWAVE)/ ;\
	cp -rv ./ecgpuwaveMakefile ./$(ECGPUWAVE)/Makefile 
ecgpuwave.isbuilt: BUILD32_DIR:=$(INSTALL_DIR)/$(BUILD_DIR)
ecgpuwave.isbuilt: wfdb32.isbuilt ecgpuwave-gfortran-addon.zip
	cd ./tmp/$(ECGPUWAVE);\
	export BUILD32_DIR;\
	make ;\
	cp -v ecgpuwave $(INSTALL_DIR)/$(BUILD_DIR)/bin ;\
	cp -v /usr/local/lib/i386/libquadmath.0.dylib $(INSTALL_DIR)/$(BUILD_DIR)/lib ;\
	cp -v /usr/local/lib/i386/libgfortran.3.dylib $(INSTALL_DIR)/$(BUILD_DIR)/lib ;\
	cp -v /usr/local/lib/i386/libgfortran.dylib $(INSTALL_DIR)/$(BUILD_DIR)/lib ;\
	cp -rv $(INSTALL_DIR)/$(BUILD_DIR32)/lib/* $(INSTALL_DIR)/$(BUILD_DIR)/lib ;\
	touch ../../ecgpuwave.isbuilt 
ecgpuwave: ecgpuwave.isbuilt


all: curl wfdb mse edr wfdb32 ecgpuwave	


#The install target will copy the contents of the build directory
#to the current lib lib64 and bin directories and remove configuration specific files
install:
	rm -rv ./bin/* ;\
	rm -rf ./lib/* ;\
	rm -rf ./lib64/* ;\
	cp -R ./$(BUILD_DIR)/bin/* ./bin/ ;\
	cp -R ./$(BUILD_DIR)/lib/* ./lib/ ;\
	cp -R ./$(BUILD_DIR)/lib64/* ./lib64/ ;\
	rm -f ./bin/psfd ./bin/hrmem ./bin/cur-config ./bin/hrfft ./bin/pschar ./bin/hrlomb ./bin/hrplot ./bin/setwfdb ./bin/cshsetwfdb ./bin/curl-config ./lib/libcurl.la ./lib/pkgconfig/libcurl.pc

#Test will run the MATLAB Unit tests headlessly 
#This requires the complete code repository and that the Makefile
#is run under [root-install-dir]/wfdb-app-tolbox/mcode/nativelibs/linux-amd64/ directory
test: 
	cd ../../ ;\
	$(MATLAB) -nodesktop -nosplash -nodisplay -r "addpath(pwd);cd ../UnitTests;BatchTest;exit"


#Package will clean all temporary files 
#generated by this Makefile so that the Toolbox binaries 
#are ready to be deployed
package: install
	rm -rf ./tmp/ ;\
	rm -rf ./$(BUILD_DIR) ;\
	rm -rf ./$(BUILD_DIR32) ;\
	rm curl* ecgpuwave* edr* mse* wfdb* *~
