#include <cstdio>
#include <cmath>

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include <arpa/inet.h>

#include <png.h>

#include "io.hpp"
#include "octree.hpp"

namespace io
{

void write_ppm(std::ostream& os, data::pixel<uint8_t>* src, uint32_t w, uint32_t h)
{
    // write the header and record that we're writing in big endian order
    os << "P6" << std::endl << w << " " << h << std::endl << 255 << std::endl;
    os.write(reinterpret_cast<char*>(src), 3*w*h);
}

namespace nbo
{

inline int32_t utos(uint32_t u)
{
	return *(reinterpret_cast<int32_t*>(&u));
}

inline int16_t utos(uint16_t u)
{
	return *(reinterpret_cast<int16_t*>(&u));
}

inline uint32_t stou(int32_t u)
{
	return *(reinterpret_cast<uint32_t*>(&u));
}

inline uint16_t stou(int16_t u)
{
	return *(reinterpret_cast<uint16_t*>(&u));
}

std::ostream& write(std::ostream& os, uint32_t val)
{
	uint32_t nval = htonl(val);
	os.write(reinterpret_cast<char*>(&nval), sizeof(uint32_t));
	return os;
}

std::ostream& write(std::ostream& os, uint16_t val)
{
	uint16_t nval = htons(val);
	os.write(reinterpret_cast<char*>(&nval), sizeof(uint16_t));
	return os;
}

std::ostream& write(std::ostream& os, int32_t val) { return write(os, stou(val)); }
std::ostream& write(std::ostream& os, int16_t val) { return write(os, stou(val)); }

std::istream& read(std::istream& is, uint32_t& r_val)
{
	uint32_t nval;
	is.read(reinterpret_cast<char*>(&nval), sizeof(uint32_t));
	r_val = ntohl(nval);
	return is;
}

std::istream& read(std::istream& is, uint16_t& r_val)
{
	uint16_t nval;
	is.read(reinterpret_cast<char*>(&nval), sizeof(uint16_t));
	r_val = ntohs(nval);
	return is;
}

std::istream& read(std::istream& is, int32_t& r_val) { uint32_t uv; read(is, uv); r_val = utos(uv); return is; }
std::istream& read(std::istream& is, int16_t& r_val) { uint16_t uv; read(is, uv); r_val = utos(uv); return is; }

}

static void abort_(const char * s, ...)
{
        va_list args;
        va_start(args, s);
        vfprintf(stderr, s, args);
        fprintf(stderr, "\n");
        va_end(args);
        abort();
}

void read_png(const char* file_name, data::image<data::pixel<float> >& im, data::image<float>& alpha)
{
	png_byte header[8];    // 8 is the maximum size that can be checked

        /* open file and test for it being a png */
        FILE *fp = fopen(file_name, "rb");
        if (!fp)
                abort_("[read_png_file] File %s could not be opened for reading", file_name);
        if(8 != fread(header, 1, 8, fp))
                abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);
        if (png_sig_cmp(header, 0, 8))
                abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);


        /* initialize stuff */
        png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[read_png_file] png_create_read_struct failed");

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[read_png_file] png_create_info_struct failed");

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during init_io");

        png_init_io(png_ptr, fp);
        png_set_sig_bytes(png_ptr, 8);

        png_read_info(png_ptr, info_ptr);

        im.width = alpha.width = png_get_image_width(png_ptr, info_ptr);
        im.height = alpha.height = png_get_image_height(png_ptr, info_ptr);

        png_byte color_type = png_get_color_type(png_ptr, info_ptr);
	if (color_type != PNG_COLOR_TYPE_RGB_ALPHA)
                abort_("[read_png_file] PNG is not RGBA");

        png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);
	if (bit_depth != 8)
                abort_("[read_png_file] PNG is not 8-bit depth");

        // int number_of_passes = png_set_interlace_handling(png_ptr);
        png_read_update_info(png_ptr, info_ptr);

        /* read file */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during read_image");

        png_bytep* row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * im.height);
        for (int y=0; y<im.height; y++)
                row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));

        png_read_image(png_ptr, row_pointers);

	// copy into images
	im.pixels.reserve(im.width * im.height);
	alpha.pixels.reserve(alpha.width * alpha.height);

        for (int y=0; y<im.height; y++)
	{
		png_byte* p_row = row_pointers[y];

		for (int x=0; x<im.width; ++x)
		{
			png_byte* p_pix = p_row + 4*x;

			float r = static_cast<float>(p_pix[0]) / 255.f;
			float g = static_cast<float>(p_pix[1]) / 255.f;
			float b = static_cast<float>(p_pix[2]) / 255.f;
			float a = static_cast<float>(p_pix[3]) / 255.f;

			im.pixels.push_back(data::pixel<float>(r*r, g*g, b*b));
			alpha.pixels.push_back(a);
		}
	}

        for (int y=0; y<im.height; y++)
		free(row_pointers[y]);
	free(row_pointers);

        fclose(fp);
}

}

std::ostream& operator << (std::ostream& os, const octree::location& loc)
{
	os << "(" << loc.x << "," << loc.y << "," << loc.z << ")";
	return os;
}
