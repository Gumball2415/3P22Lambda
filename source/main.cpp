#include <filesystem>

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

		cmd.add(inpath_arg);
		cmd.add(outpath_arg);

		cmd.parse(argc, argv);

		inpath = inpath_arg.getValue();
		outpath = outpath_arg.getValue();

		// check
		if (inpath.is_absolute())

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

	// step 2: process each frame through ProcessFrame()

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
	PhosphorFrame P22(width, height, bitdepth, 0.5, 0.5, 0.5);

	for (const auto &frame : framepaths)
		P22.ProcessFrame(frame.string(), outpath.string());
	// step 3: put rendered frames into output folder

	// step 4: profit
	return 0;
}

PhosphorFrame::PhosphorFrame(uint32_t frame_width, uint32_t frame_height, uint32_t frame_bitdepth, double rp22_falloff, double gp22_falloff, double bp22_falloff)
{
	Width = frame_width;
	Height = frame_height;
	Bitdepth = frame_bitdepth;
	RP22Falloff = rp22_falloff;
	GP22Falloff = gp22_falloff;
	BP22Falloff = bp22_falloff;
	RP22FrameBuffer.resize(frame_width * frame_height, 0);
	GP22FrameBuffer.resize(frame_width * frame_height, 0);
	BP22FrameBuffer.resize(frame_width * frame_height, 0);
	RP22FrameBufferPrev.resize(frame_width * frame_height, 0);
	GP22FrameBufferPrev.resize(frame_width * frame_height, 0);
	BP22FrameBufferPrev.resize(frame_width * frame_height, 0);
	FrameBufferFinal.resize(frame_width * frame_height, 0);
}

PhosphorFrame::~PhosphorFrame()
{
}

void PhosphorFrame::ProcessFrame(std::string file_path, std::string output_directory)
{
	std::filesystem::path outputfile = std::filesystem::path(output_directory) / std::filesystem::path(file_path).filename();
	std::cout << outputfile << std::endl;

	// decode input .png into raw framebuffer
	std::vector<uint8_t> rawimage;

	// convert raw image into framebuffer
	// todo: spawn new threads to speed this up
	for (int i = 0; i < rawimage.size(); i++) {
		// todo: convert raw image RGB to SMPTE colorimetry stuff, sRGB assumed
		
		RP22FrameBuffer.at(i) = ProcessPixel((static_cast<double>(rawimage.at(i)) / static_cast<double>(Bitdepth)), RP22Falloff, RP22FrameBufferPrev.at(i));
		GP22FrameBuffer.at(i) = ProcessPixel((static_cast<double>(rawimage.at(i)) / static_cast<double>(Bitdepth)), GP22Falloff, GP22FrameBufferPrev.at(i));
		BP22FrameBuffer.at(i) = ProcessPixel((static_cast<double>(rawimage.at(i)) / static_cast<double>(Bitdepth)), BP22Falloff, BP22FrameBufferPrev.at(i));

	}

	for (int i = 0; i < FrameBufferFinal.size(); i++) {
		// todo: do colorimetry stuff to convert back to sRGB
	}

	RP22FrameBufferPrev = RP22FrameBuffer;
	GP22FrameBufferPrev = GP22FrameBuffer;
	BP22FrameBufferPrev = BP22FrameBuffer;
}

double PhosphorFrame::ProcessPixel(double input_pixel, double phosphor_falloff, double previous_pixel)
{
	return std::max(previous_pixel * phosphor_falloff, input_pixel);
}
