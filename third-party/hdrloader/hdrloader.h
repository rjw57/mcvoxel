#ifndef MC_VOXEL_THIRD_PARTY_HDRLOADER_HDRLOADER_H__
#define MC_VOXEL_THIRD_PARTY_HDRLOADER_HDRLOADER_H__

/***********************************************************************************
	Created:	17:9:2002
	FileName: 	hdrloader.h
	Author:		Igor Kravtchenko
	
	Info:		Load HDR image and convert to a set of float32 RGB triplet.
************************************************************************************/

class HDRLoaderResult {
public:
	int width, height;
	// each pixel takes 3 float32, each component can be of any value...
	float *cols;
};

class HDRLoader {
public:
	static bool load(const char *fileName, HDRLoaderResult &res);
};

#endif // MC_VOXEL_THIRD_PARTY_HDRLOADER_HDRLOADER_H__
