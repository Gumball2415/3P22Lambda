#include <filesystem>
#include <sstream>

#include "main.h"
#include "config.h"
#include "lodepng/lodepng.h"
#include <tclap/CmdLine.h>

int main(int argc, char* argv[])
{
	std::cout << PROJECT_NAME << " v." << PROJECT_VERSION << std::endl
		<< PROJECT_DESCRIPTION << std::endl
		<< "uses LodePNG, Copyright (c) 2005-2010 Lode Vandevenne" << std::endl;
	
	std::filesystem::path inpath, outpath;
	std::string fpsstring;

	try {
		TCLAP::CmdLine cmd("uses TCLAP, Copyright (c) 2017-2021 Google LLC, (c) 2012 - 2016 Daniel Aarno, (c) 2003 - 2012 Michael E.Smoot", ' ', PROJECT_VERSION);
		TCLAP::ValueArg<std::string> inpath_arg(
			"i",
			"input",
			"input path to folder of .png sequence, preferrably sorted",
			true,
			"",
			"path");
		TCLAP::ValueArg<std::string> outpath_arg(
			"o",
			"output",
			"output path to folder of .png sequence",
			true,
			"",
			"path");
		TCLAP::ValueArg<std::string> fps_arg(
			"f",
			"fps",
			"frames per second of .png sequence",
			true,
			"60/1",
			"framerate");

		cmd.add(inpath_arg);
		cmd.add(outpath_arg);
		cmd.add(fps_arg);

		cmd.parse(argc, argv);

		inpath = inpath_arg.getValue();
		outpath = outpath_arg.getValue();
		fpsstring = fps_arg.getValue();

		std::cout << "input folder: " << inpath << std::endl;
		std::cout << "output folder: " << outpath << std::endl;
		
	}
	catch (TCLAP::ArgException& e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

	// step 1: count input .png frames
	std::vector<std::filesystem::path> framepaths;
	for (const auto& p : std::filesystem::directory_iterator(inpath)) {
		if (p.path().extension() == ".png")
			framepaths.push_back(p);
	}
	std::cout << "frame count: " << framepaths.size() << std::endl;

	// step 2 + 3: process each frame through ProcessFrame()
	// and put rendered frames into output folder

	// determine height and width by decoding the first frame
	std::string previewframe_filepath = framepaths.at(0).string();
	std::vector<uint8_t> previewframe_buffer;
	std::vector<uint8_t> previewframe_image;
	lodepng::State previewframe_state;
	uint32_t width = 0, height = 0;
	lodepng::load_file(previewframe_buffer, previewframe_filepath);

	unsigned int error = lodepng::decode(previewframe_image, width, height, previewframe_state, previewframe_buffer);
	if (error) {
		std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
		return 0;
	}
	uint32_t bitdepth = previewframe_state.info_png.color.bitdepth;

	std::cout << "width: " << width << std::endl
		<< "height: " << height << std::endl
		<< "bitdepth: " << bitdepth << std::endl;

	double fps, fpsnum, fpsden;
	char c;
	std::stringstream fpsstream(fpsstring);
	fpsstream >> fpsnum;
	if (fpsstream.fail() || fpsstream.get(c)) {
		// if "/" parsed. then it's a fraction
		fpsstream >> fpsden;
		fps = fpsnum / fpsden;
	}
	else {
		fps = fpsnum;
	}


	PhosphorFrame P22(width, height, fps);

	std::cout << std::fixed << std::setprecision(3);
	for (int i = 0; i < framepaths.size(); i++) {
		std::cout << "\rprogress: " << i + 1 << " / " << framepaths.size()
			<< " (" << static_cast<double>(i + 1.0) / static_cast<double>(framepaths.size()) * 100.0 << "%)" << std::flush;
		P22.ProcessFrame(framepaths.at(i).string(), outpath.string());
	}

	// step 4: profit
	std::cout << std::endl << "done ^w^" << std::endl;
	return 0;
}

PhosphorFrame::PhosphorFrame(uint32_t frame_width, uint32_t frame_height, double fps_delta)
{
	Width = frame_width;
	Height = frame_height;

	// derived from luminance components of sRGB
	RP22Falloff = (0.2126 / 4.0) + 0.4;
	GP22Falloff = (0.7152 / 4.0) + 0.4;
	BP22Falloff = (0.0722 / 4.0) + 0.4;
	FPS = fps_delta;
	P22ShortFrameBuffer.resize(frame_width * frame_height, std::vector<double>(3, 0.0));
	P22LongFrameBuffer.resize(frame_width * frame_height, std::vector<double>(3, 0.0));
}

PhosphorFrame::~PhosphorFrame()
{
}

void PhosphorFrame::ProcessFrame(std::string file_path, std::string output_directory)
{
	// decode input .png into raw framebuffer

	std::vector<uint8_t> rawimage;
	std::vector<uint8_t> buffer;
	unsigned w, h;
	lodepng::State state;
	state.info_raw.bitdepth = 16; // higher bit depth for more accurate color math
	state.info_raw.colortype = LCT_RGBA; // in case we are dealing with alpha

	lodepng::load_file(buffer, file_path);
	unsigned error = lodepng::decode(rawimage, w, h, state, buffer);

	if (error) {
		std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
		return;
	}
	if (w != Width || h != Height) {
		std::cout << "error: " << "at " << file_path << std::endl
			<< "input image dimensions doesn't match first frame!" << std::endl;
		return;
	}

	FrameBufferFinal.clear();
	FrameBufferFinal.resize(rawimage.size());

	// convert raw image into framebuffer
	for (int i = 0; i < rawimage.size(); i += 8) {
		// todo: spawn new threads to speed this up
		{
			// raw image format is 0xRRRR 0xGGGG 0xBBBB 0xAAAA, big endian
			uint16_t pixel_r = (rawimage[i + 0] << 8) | rawimage[i + 1];
			uint16_t pixel_g = (rawimage[i + 2] << 8) | rawimage[i + 3];
			uint16_t pixel_b = (rawimage[i + 4] << 8) | rawimage[i + 5];

			// convert raw image RGB to CIE XYZ coordinates of P22 primaries, sRGB assumed
			std::vector<double> pixel_convert_p22 = sRGBToP22(pixel_r, pixel_g, pixel_b);

			// process pixel
			P22ShortFrameBuffer[i / 8] = ProcessPixel(pixel_convert_p22, P22ShortFrameBuffer[i / 8], RP22Falloff, GP22Falloff, BP22Falloff);
			P22LongFrameBuffer[i / 8] = ProcessPixel(pixel_convert_p22, P22LongFrameBuffer[i / 8], 0.9, 0.9, 0.9);

			// convert back to sRGB
			std::vector<uint16_t> pixel_convert_srgb = P22TosRGB(
				(P22LongFrameBuffer[i / 8][0] * 0.05) + P22ShortFrameBuffer[i / 8][0],
				(P22LongFrameBuffer[i / 8][1] * 0.05) + P22ShortFrameBuffer[i / 8][1],
				(P22LongFrameBuffer[i / 8][2] * 0.05) + P22ShortFrameBuffer[i / 8][2]);

			FrameBufferFinal[i + 0] = (pixel_convert_srgb[0] & 0xFF00) >> 8;
			FrameBufferFinal[i + 1] = (pixel_convert_srgb[0] & 0xFF);
			FrameBufferFinal[i + 2] = (pixel_convert_srgb[1] & 0xFF00) >> 8;
			FrameBufferFinal[i + 3] = (pixel_convert_srgb[1] & 0xFF);
			FrameBufferFinal[i + 4] = (pixel_convert_srgb[2] & 0xFF00) >> 8;
			FrameBufferFinal[i + 5] = (pixel_convert_srgb[2] & 0xFF);
			// preserve transparency
			FrameBufferFinal[i + 6] = rawimage[i + 6];
			FrameBufferFinal[i + 7] = rawimage[i + 7];
		}
	}

	// save frame
	std::filesystem::path outputfile = std::filesystem::path(output_directory) / std::filesystem::path(file_path).filename();

	// use the same settings as the original
	buffer.clear();
	// new colors are generated using this filter, so force RGB/RGBA
	if (state.info_png.color.colortype != LCT_RGB || state.info_png.color.colortype != LCT_RGBA) {
		if (state.info_png.color.colortype == LCT_GREY_ALPHA)
			state.info_png.color.colortype = LCT_RGBA;
		else
			state.info_png.color.colortype = LCT_RGB;
	}

	error = lodepng::encode(buffer, FrameBufferFinal, w, h, state);

	if (error) {
		std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
		return;
	}
	else
		lodepng::save_file(buffer, outputfile.string());
}

std::vector<double> PhosphorFrame::ProcessPixel(std::vector<double> in, std::vector<double> prev, double rp22_falloff, double gp22_falloff, double bp22_falloff)
{

	double after_r = (in[0] + prev[0] * 0.95);
	double after_g = (in[1] + prev[1] * 0.95);
	double after_b = (in[2] + prev[2] * 0.95);

	return {
		(in[0] * (1.0 - rp22_falloff)) + (prev[0] * rp22_falloff),
		(in[1] * (1.0 - gp22_falloff)) + (prev[1] * gp22_falloff),
		(in[2] * (1.0 - bp22_falloff)) + (prev[2] * bp22_falloff)
	};
}

std::vector<double> PhosphorFrame::sRGBToP22(uint16_t pixel_r, uint16_t pixel_g, uint16_t pixel_b)
{
	double rgb[3] = {
		static_cast<double>(pixel_r) / static_cast<double>(0xFFFF),
		static_cast<double>(pixel_g) / static_cast<double>(0xFFFF),
		static_cast<double>(pixel_b) / static_cast<double>(0xFFFF)
	};

	// sRGB to linear RGB
	// reference: https://www.color.org/chardata/rgb/srgb.xalter
	for (int i = 0; i < 3; i++) {
		rgb[i] = (rgb[i] <= 0.04045) ? rgb[i] / 12.92 : std::pow(((rgb[i] + 0.055) / 1.055), 2.4);
	}

	return { rgb[0], rgb[1], rgb[2] };
}

std::vector<uint16_t> PhosphorFrame::P22TosRGB(double pixel_r, double pixel_g, double pixel_b)
{
	double rgb[3] = { pixel_r, pixel_g, pixel_b };

	// optional: interpret linear RGB as NTSC 1963 primaries

	// linear RGB to sRGB
	// reference: https://www.color.org/chardata/rgb/srgb.xalter
	for (int i = 0; i < 3; i++) {
		rgb[i] = (rgb[i] <= 0.0031308) ? rgb[i] * 12.92 : (1.055 * std::pow(rgb[i], (1.0 / 2.4))) - 0.055;

		// clamp RGB
		if (rgb[i] <= 0.0) rgb[i] = 0.0;
		if (rgb[i] >= 1.0) rgb[i] = 1.0;
	}

	return {
		static_cast<uint16_t>(rgb[0] * static_cast<double>(0xFFFF)),
		static_cast<uint16_t>(rgb[1] * static_cast<double>(0xFFFF)),
		static_cast<uint16_t>(rgb[2] * static_cast<double>(0xFFFF))
	};
}
