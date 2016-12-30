#
#
#
noiseModelApp: noiseModelApp.C
	g++ -std=c++0x -L$(ANITA_LIB_DIR) noiseModelApp.C -o noiseModelApp `root-config --cflags --glibs` -I$(ANITA_UTIL_INSTALL_DIR)/include -lRootFftwWrapper

.PHONY: clean

clean: 
	rm -f noiseModelApp
