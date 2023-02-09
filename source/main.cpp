#include <filesystem>
#include <sstream>
#include <thread>

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
	double fps;
	double rp22_short, gp22_short, bp22_short;
	double rp22_long, gp22_long, bp22_long;
	double p22_long_mix;
	bool ntsc_1953, multithreaded;

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

		TCLAP::SwitchArg old_primaries_arg(
			"p",
			"primaries",
			"use NTSC 1953 primaries",
			false
		);
		TCLAP::SwitchArg multithreaded_arg(
			"m",
			"multithread",
			"use multithreading",
			false
		);

		// derived from luminance components of sRGB
		TCLAP::ValueArg<double> r_short_falloff_arg(
			"1",
			"rshort",
			"red phosphor short falloff rate",
			false,
			(0.2126 / 4.0),
			"0.0 - 1.0 float value");
		TCLAP::ValueArg<double> g_short_falloff_arg(
			"2",
			"gshort",
			"green phosphor short falloff rate",
			false,
			(0.7152 / 4.0),
			"0.0 - 1.0 float value");
		TCLAP::ValueArg<double> b_short_falloff_arg(
			"3",
			"bshort",
			"blue phosphor short falloff rate",
			false,
			(0.0722 / 4.0),
			"0.0 - 1.0 float value");

		TCLAP::ValueArg<double> r_long_falloff_arg(
			"4",
			"rlong",
			"red phosphor long falloff rate",
			false,
			(0.2126 / 10.0) + 0.1,
			"0.0 - 1.0 float value");
		TCLAP::ValueArg<double> g_long_falloff_arg(
			"5",
			"glong",
			"green phosphor long falloff rate",
			false,
			(0.7152 / 10.0) + 0.8,
			"0.0 - 1.0 float value");
		TCLAP::ValueArg<double> b_long_falloff_arg(
			"6",
			"blong",
			"blue phosphor long falloff rate",
			false,
			(0.0722 / 10.0) + 0.9,
			"0.0 - 1.0 float value");

		TCLAP::ValueArg<double> long_falloff_mix(
			"7",
			"longmix",
			"long phosphor falloff mix",
			false,
			0.01,
			"0.0 - 1.0 float value");

		cmd.add(inpath_arg);
		cmd.add(outpath_arg);
		cmd.add(fps_arg);
		cmd.add(old_primaries_arg);
		cmd.add(multithreaded_arg);

		cmd.add(r_short_falloff_arg);
		cmd.add(g_short_falloff_arg);
		cmd.add(b_short_falloff_arg);
		cmd.add(r_long_falloff_arg);
		cmd.add(g_long_falloff_arg);
		cmd.add(b_long_falloff_arg);
		cmd.add(long_falloff_mix);

		cmd.parse(argc, argv);

		inpath = inpath_arg.getValue();
		outpath = outpath_arg.getValue();
		fpsstring = fps_arg.getValue();
		ntsc_1953 = old_primaries_arg.getValue();
		multithreaded = multithreaded_arg.getValue();

		rp22_short = r_short_falloff_arg.getValue();
		gp22_short = g_short_falloff_arg.getValue();
		bp22_short = b_short_falloff_arg.getValue();

		rp22_long = r_long_falloff_arg.getValue();
		gp22_long = g_long_falloff_arg.getValue();
		bp22_long = b_long_falloff_arg.getValue();
		p22_long_mix = long_falloff_mix.getValue();

		char c;
		double fpsnum, fpsden;
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

		std::cout << std::fixed << std::setprecision(4);
		std::cout << "stats:" << std::endl
			<< "    input folder:               " << inpath << std::endl
			<< "    output folder:              " << outpath << std::endl
			<< "    calculated FPS:             " << fps << std::endl
			<< "    use NTSC 1953 primaries:    " << (ntsc_1953 ? "true" : "false") << std::endl
			<< "    use multithreading:         " << (multithreaded ? "true" : "false") << std::endl << std::endl
			<< "    red P22 short falloff:      " << rp22_short << std::endl
			<< "    green P22 short falloff:    " << gp22_short << std::endl
			<< "    blue P22 short falloff:     " << bp22_short << std::endl << std::endl
			<< "    red P22 long falloff:       " << rp22_long << std::endl
			<< "    green P22 long falloff:     " << gp22_long << std::endl
			<< "    blue P22 long falloff:      " << bp22_long << std::endl
			<< "    P22 long falloff mix:       " << p22_long_mix << std::endl << std::endl;
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

	if (framepaths.size() <= 0) {
		std::cout << "error: no .png frames detected in input folder" << std::endl;
		return 0;
	}


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

	std::cout << "frame stats:" << std::endl
		<< "    frame count:            " << framepaths.size() << std::endl
		<< "    width:                  " << width << std::endl
		<< "    height:                 " << height << std::endl
		<< "    bitdepth per channel:   " << bitdepth << std::endl << std::endl;


	PhosphorFrame P22(width, height, fps, rp22_short, gp22_short, bp22_short, rp22_long, gp22_long, bp22_long, p22_long_mix, ntsc_1953, multithreaded);

	for (int i = 0; i < framepaths.size(); i++) {
		std::cout << "\rprogress: " << i + 1 << " / " << framepaths.size()
			<< " (" << static_cast<double>(i + 1.0) / static_cast<double>(framepaths.size()) * 100.0 << "%)" << std::flush;
		P22.ProcessFrame(framepaths.at(i).string(), outpath.string());
	}

	// step 4: profit
	std::cout << std::endl << "done ^w^" << std::endl;
	return 0;
}

PhosphorFrame::PhosphorFrame(uint32_t frame_width, uint32_t frame_height, double fps_delta, double rp22_short, double gp22_short, double bp22_short, double rp22_long, double gp22_long, double bp22_long, double p22_long_mix, bool old_primaries, bool multithreading)
{
	Width = frame_width;
	Height = frame_height;

	RP22ShortFalloff = rp22_short;
	GP22ShortFalloff = gp22_short;
	BP22ShortFalloff = bp22_short;

	RP22LongFalloff = rp22_long;
	GP22LongFalloff = gp22_long;
	BP22LongFalloff = bp22_long;

	P22LongFalloffMix = p22_long_mix;
	
	FPS = fps_delta;

	UseNTSC1953Primaries = old_primaries;
	UseMultithreading = multithreading;

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

	unsigned int threadcount = std::thread::hardware_concurrency();

	if (threadcount != 0 && UseMultithreading) {
		unsigned int i = 0;

		unsigned int framechunk, framechunk_r = 0;
		framechunk = (Width * Height) / threadcount;

		if (((Width * Height) % threadcount))
			framechunk_r = (Width * Height) - (framechunk * threadcount);
		
		// spawn new threads to speed this up
		std::vector<std::thread> threadvectorpixel;

		for (int j = 0; j < threadcount; j++) {
			if (framechunk_r && (j == (threadcount - 1))) {
				threadvectorpixel.emplace_back([&, i]() { ThreadProcessPixels(rawimage, i, i + framechunk_r - 1); });
				i += framechunk_r;
				break;
			}
			threadvectorpixel.emplace_back([&, i]() { ThreadProcessPixels(rawimage, i, i + framechunk - 1); });
			i += framechunk;
		}
		_ASSERT((i * 8) == rawimage.size());
		for (auto& thread : threadvectorpixel)
			thread.join();
	}
	else {
		ThreadProcessPixels(std::ref(rawimage), 0, rawimage.size());
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

void PhosphorFrame::ThreadProcessPixels(std::vector<uint8_t> &rawimage_buffer, unsigned int start_pixel, unsigned int end_pixel)
{
	for (int i = start_pixel; i <= end_pixel; i++) {
		unsigned int buffer_index = (i * 8);
		// raw image format is 0xRRRR 0xGGGG 0xBBBB 0xAAAA, big endian
		uint16_t pixel_r = (rawimage_buffer[buffer_index + 0] << 8) | rawimage_buffer[buffer_index + 1];
		uint16_t pixel_g = (rawimage_buffer[buffer_index + 2] << 8) | rawimage_buffer[buffer_index + 3];
		uint16_t pixel_b = (rawimage_buffer[buffer_index + 4] << 8) | rawimage_buffer[buffer_index + 5];

		// convert raw image RGB to CIE XYZ coordinates of P22 primaries, sRGB assumed
		std::vector<double> pixel_convert_p22 = sRGBToP22(pixel_r, pixel_g, pixel_b);

		// process pixel
		P22ShortFrameBuffer[i] = ProcessPixel(pixel_convert_p22, P22ShortFrameBuffer[i], RP22ShortFalloff, GP22ShortFalloff, BP22ShortFalloff, FPS);
		P22LongFrameBuffer[i] = ProcessPixel(pixel_convert_p22, P22LongFrameBuffer[i], RP22LongFalloff, GP22LongFalloff, BP22LongFalloff, FPS);

		double result[3] = {
			std::max((P22LongFrameBuffer[i][0] * P22LongFalloffMix), P22ShortFrameBuffer[i][0]),
			std::max((P22LongFrameBuffer[i][1] * P22LongFalloffMix), P22ShortFrameBuffer[i][1]),
			std::max((P22LongFrameBuffer[i][2] * P22LongFalloffMix), P22ShortFrameBuffer[i][2]),
		};

		// convert back to sRGB
		std::vector<uint16_t> pixel_convert_srgb = P22TosRGB(result[0], result[1], result[2]);

		FrameBufferFinal[buffer_index + 0] = (pixel_convert_srgb[0] & 0xFF00) >> 8;
		FrameBufferFinal[buffer_index + 1] = (pixel_convert_srgb[0] & 0xFF);
		FrameBufferFinal[buffer_index + 2] = (pixel_convert_srgb[1] & 0xFF00) >> 8;
		FrameBufferFinal[buffer_index + 3] = (pixel_convert_srgb[1] & 0xFF);
		FrameBufferFinal[buffer_index + 4] = (pixel_convert_srgb[2] & 0xFF00) >> 8;
		FrameBufferFinal[buffer_index + 5] = (pixel_convert_srgb[2] & 0xFF);
		// preserve transparency
		FrameBufferFinal[buffer_index + 6] = rawimage_buffer[buffer_index + 6];
		FrameBufferFinal[buffer_index + 7] = rawimage_buffer[buffer_index + 7];
	}
}

std::vector<double> PhosphorFrame::ProcessPixel(std::vector<double> &in, std::vector<double> &prev, double &rp22_falloff, double &gp22_falloff, double &bp22_falloff, double &fps)
{
	// todo: simulate phosphor afterglow entanglement
	std::vector<double> output = {
		// default falloff values were optimized for 60fps video
		std::max(in[0], (prev[0] * (rp22_falloff / 60.0 * fps))),
		std::max(in[1], (prev[1] * (gp22_falloff / 60.0 * fps))),
		std::max(in[2], (prev[2] * (bp22_falloff / 60.0 * fps))),
	};
	return output;
}

std::vector<double> PhosphorFrame::sRGBToP22(uint16_t &pixel_r, uint16_t &pixel_g, uint16_t &pixel_b)
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

std::vector<uint16_t> PhosphorFrame::P22TosRGB(double &pixel_r, double &pixel_g, double &pixel_b)
{
	double rgb[3] = { pixel_r, pixel_g, pixel_b };

	// optional: interpret linear RGB as NTSC 1963 primaries
	// thanks NewRisingSun!
	if (UseNTSC1953Primaries) {
		// ignoring whitepoint difference
		double ntsc_1953_rgb[3] = {
			 1.4607 * rgb[0] - 0.3845 * rgb[1] - 0.0761 * rgb[2],
			-0.0266 * rgb[0] + 0.9654 * rgb[1] + 0.0612 * rgb[2],
			-0.0264 * rgb[0] + 0.0414 * rgb[1] + 1.0678 * rgb[2]
		};
		rgb[0] = ntsc_1953_rgb[0];
		rgb[1] = ntsc_1953_rgb[1];
		rgb[2] = ntsc_1953_rgb[2];
	}

	// linear RGB to sRGB
	// reference: https://www.color.org/chardata/rgb/srgb.xalter
	for (int i = 0; i < 3; i++) {
		rgb[i] = (rgb[i] <= 0.0031308) ? rgb[i] * 12.92 : (1.055 * std::pow(rgb[i], (1.0 / 2.4))) - 0.055;

		// clamp RGB
		if (rgb[i] <= 0.0) rgb[i] = 0.0;
		if (rgb[i] >= 1.0) rgb[i] = 1.0;
	}

	return {
		static_cast<uint16_t>(rgb[0] * 0xFFFF),
		static_cast<uint16_t>(rgb[1] * 0xFFFF),
		static_cast<uint16_t>(rgb[2] * 0xFFFF)
	};
}
