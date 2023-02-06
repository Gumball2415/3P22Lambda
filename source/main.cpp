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
	double fps;
	bool srgb_direct_transfer, p22_direct_transfer;

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
		TCLAP::ValueArg<double> fps_arg(
			"f",
			"fps",
			"frames per second of .png sequence",
			true,
			0.0,
			"framerate");
		TCLAP::SwitchArg srgb_direct_arg(
			"s",
			"srgbdirect",
			"directly transfer RGB to P22 gamut coordinates",
			false
		);
		TCLAP::SwitchArg p22_direct_arg(
			"p",
			"p22direct",
			"directly transfer P22 to sRGB gamut coordinates",
			false
		);

		cmd.add(inpath_arg);
		cmd.add(outpath_arg);
		cmd.add(fps_arg);
		cmd.add(srgb_direct_arg);
		cmd.add(p22_direct_arg);

		cmd.parse(argc, argv);

		inpath = inpath_arg.getValue();
		outpath = outpath_arg.getValue();
		fps = fps_arg.getValue();
		srgb_direct_transfer = srgb_direct_arg.getValue();
		p22_direct_transfer = p22_direct_arg.getValue();

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
	PhosphorFrame P22(width, height, 0.5, 0.4, 0.3, fps);

	std::cout << std::fixed << std::setprecision(3);
	for (int i = 0; i < framepaths.size(); i++) {
		std::cout << "\rprogress: " << i + 1 << " / " << framepaths.size()
			<< " (" << static_cast<double>(i + 1.0) / static_cast<double>(framepaths.size()) * 100.0 << "%)" << std::flush;
		P22.ProcessFrame(framepaths.at(i).string(), outpath.string(), srgb_direct_transfer, p22_direct_transfer);
	}

	// step 4: profit
	std::cout << std::endl << "done ^w^" << std::endl;
	return 0;
}

PhosphorFrame::PhosphorFrame(uint32_t frame_width, uint32_t frame_height, double rp22_falloff, double gp22_falloff, double bp22_falloff, double fps_delta)
{
	Width = frame_width;
	Height = frame_height;
	RP22Falloff = rp22_falloff;
	GP22Falloff = gp22_falloff;
	BP22Falloff = bp22_falloff;
	FPSDelta = fps_delta;
	RP22FrameBuffer.resize(frame_width * frame_height, { 0.3127, 0.3290, 0.0 });
	GP22FrameBuffer.resize(frame_width * frame_height, { 0.3127, 0.3290, 0.0 });
	BP22FrameBuffer.resize(frame_width * frame_height, { 0.3127, 0.3290, 0.0 });
}

PhosphorFrame::~PhosphorFrame()
{
}

void PhosphorFrame::ProcessFrame(std::string file_path, std::string output_directory, bool srgb_direct_transfer, bool p22_direct_transfer)
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

	FrameBufferFinal.resize(rawimage.size());

	// convert raw image into framebuffer
	for (int i = 0; i < rawimage.size(); i += 8) {
		// todo: spawn new threads to speed this up
		{
			// raw image format is 0xRRRR 0xGGGG 0xBBBB 0xAAAA, big endian
			uint16_t pixel_r = (rawimage.at(i + 0) << 8) | rawimage.at(i + 1);
			uint16_t pixel_g = (rawimage.at(i + 2) << 8) | rawimage.at(i + 3);
			uint16_t pixel_b = (rawimage.at(i + 4) << 8) | rawimage.at(i + 5);

			// convert raw image RGB to CIE XYZ coordinates of P22 primaries, sRGB assumed
			std::vector<std::vector<double>> pixel_convert_p22 = sRGBToP22(pixel_r, pixel_g, pixel_b, srgb_direct_transfer);

			// process pixel
			// RP22 xyY
			RP22FrameBuffer.at(i / 8) = ProcessPixel(pixel_convert_p22.at(0), RP22Falloff, RP22FrameBuffer.at(i / 8));

			// GP22 xyY
			GP22FrameBuffer.at(i / 8) = ProcessPixel(pixel_convert_p22.at(1), GP22Falloff, GP22FrameBuffer.at(i / 8));

			// BP22 xyY
			BP22FrameBuffer.at(i / 8) = ProcessPixel(pixel_convert_p22.at(2), BP22Falloff, BP22FrameBuffer.at(i / 8));

			// convert back to sRGB
			std::vector<uint16_t> pixel_convert_srgb = P22TosRGB(RP22FrameBuffer.at(i / 8), GP22FrameBuffer.at(i / 8), BP22FrameBuffer.at(i / 8), p22_direct_transfer);

			FrameBufferFinal.at(i + 0) = (pixel_convert_srgb.at(0) & 0xFF00) >> 8;
			FrameBufferFinal.at(i + 1) = (pixel_convert_srgb.at(0) & 0xFF);
			FrameBufferFinal.at(i + 2) = (pixel_convert_srgb.at(1) & 0xFF00) >> 8;
			FrameBufferFinal.at(i + 3) = (pixel_convert_srgb.at(1) & 0xFF);
			FrameBufferFinal.at(i + 4) = (pixel_convert_srgb.at(2) & 0xFF00) >> 8;
			FrameBufferFinal.at(i + 5) = (pixel_convert_srgb.at(2) & 0xFF);
			// preserve transparency
			FrameBufferFinal.at(i + 6) = rawimage.at(i + 6);
			FrameBufferFinal.at(i + 7) = rawimage.at(i + 7);
		}
	}

	// save frame
	std::filesystem::path outputfile = std::filesystem::path(output_directory) / std::filesystem::path(file_path).filename();

	// use the same settings as the original
	buffer.clear();
	// new colors are generated using this filter, so force RGB/RGBA
	if (state.info_png.color.colortype != LCT_RGB || state.info_png.color.colortype != LCT_RGBA)
		if (state.info_png.color.colortype == LCT_GREY_ALPHA)
			state.info_png.color.colortype = LCT_RGBA;
		else
			state.info_png.color.colortype = LCT_RGB;

	error = lodepng::encode(buffer, FrameBufferFinal, w, h, state);

	if (error) {
		std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
		return;
	}
	else
		lodepng::save_file(buffer, outputfile.string());
}

std::vector<double> PhosphorFrame::ProcessPixel(std::vector<double> in, double phosphor_falloff, std::vector<double> prev)
{

	// out = in * falloff + prev * falloff
	// https://en.wikipedia.org/wiki/CIE_1931_color_space#Mixing_colors_specified_with_the_CIE_xy_chromaticity_diagram
	double out[3] = {
		((((in[0] / in[1]) * in[2]) + ((prev[0] / prev[1]) * prev[2])) / ((in[2] / in[1]) + (prev[2] / prev[1]))),
		((in[2] + prev[2]) / ((in[2] / in[1]) + (prev[2] / prev[1]))),
		(in[2] * (1.0 - phosphor_falloff) + (prev[2] * phosphor_falloff))
	};

	// usually when xy or Y is 0
	if (std::isnan(out[0])) out[0] = 0.3127;
	if (std::isnan(out[1])) out[1] = 0.3290;

	return { out[0], out[1], out[2] };
}

std::vector<std::vector<double>> PhosphorFrame::sRGBToP22(uint16_t pixel_r, uint16_t pixel_g, uint16_t pixel_b, bool direct_transfer)
{
	double rgb[3] = { static_cast<double>(pixel_r) / static_cast<double>(0xFFFF),
		static_cast<double>(pixel_g) / static_cast<double>(0xFFFF),
		static_cast<double>(pixel_b) / static_cast<double>(0xFFFF)
	};

	// sRGB to CIE XYZ to CIE xyY
	// https://en.wikipedia.org/wiki/SRGB#From_sRGB_to_CIE_XYZ

	// sRGB to linear RGB
	for (int i = 0; i < 3; i++) {
		if (rgb[i] <= 0.04045) rgb[i] /= 12.92;
		else rgb[i] = std::pow(((rgb[i] + 0.055) / 1.055), 2.4);
	}

	// linear RGB to CIE XYZ via matrix multiplication
	double srgb2cieXYZ[3][3] = {
		{0.4124, 0.3576, 0.1805},
		{0.2126, 0.7152, 0.0722},
		{0.0193, 0.1192, 0.9505}
	};

	double p22[3][3];
	if (pixel_r == 0 && pixel_g == 0 && pixel_b == 0) {
		// RP22 XYZ
		p22[0][0] = srgb2cieXYZ[0][0] + srgb2cieXYZ[0][1] + srgb2cieXYZ[0][2];
		p22[0][1] = rgb[0] * srgb2cieXYZ[1][0] + rgb[0] * srgb2cieXYZ[1][1] + rgb[0] * srgb2cieXYZ[1][2];
		p22[0][2] = srgb2cieXYZ[2][0] + srgb2cieXYZ[2][1] + srgb2cieXYZ[2][2];
		// GP22 XYZ
		p22[1][0] = srgb2cieXYZ[0][0] + srgb2cieXYZ[0][1] + srgb2cieXYZ[0][2];
		p22[1][1] = rgb[1] * srgb2cieXYZ[1][0] + rgb[1] * srgb2cieXYZ[1][1] + rgb[1] * srgb2cieXYZ[1][2];
		p22[1][2] = srgb2cieXYZ[2][0] + srgb2cieXYZ[2][1] + srgb2cieXYZ[2][2];
		// BP22 XYZ
		p22[2][0] = srgb2cieXYZ[0][0] + srgb2cieXYZ[0][1] + srgb2cieXYZ[0][2];
		p22[2][1] = rgb[2] * srgb2cieXYZ[1][0] + rgb[2] * srgb2cieXYZ[1][1] + rgb[2] * srgb2cieXYZ[1][2];
		p22[2][2] = srgb2cieXYZ[2][0] + srgb2cieXYZ[2][1] + srgb2cieXYZ[2][2];
	}
	else {
		// RP22 XYZ
		p22[0][0] = rgb[0] * srgb2cieXYZ[0][0] + rgb[0] * srgb2cieXYZ[0][1] + rgb[0] * srgb2cieXYZ[0][2];
		p22[0][1] = rgb[0] * srgb2cieXYZ[1][0] + rgb[0] * srgb2cieXYZ[1][1] + rgb[0] * srgb2cieXYZ[1][2];
		p22[0][2] = rgb[0] * srgb2cieXYZ[2][0] + rgb[0] * srgb2cieXYZ[2][1] + rgb[0] * srgb2cieXYZ[2][2];
		// GP22 XYZ
		p22[1][0] = rgb[1] * srgb2cieXYZ[0][0] + rgb[1] * srgb2cieXYZ[0][1] + rgb[1] * srgb2cieXYZ[0][2];
		p22[1][1] = rgb[1] * srgb2cieXYZ[1][0] + rgb[1] * srgb2cieXYZ[1][1] + rgb[1] * srgb2cieXYZ[1][2];
		p22[1][2] = rgb[1] * srgb2cieXYZ[2][0] + rgb[1] * srgb2cieXYZ[2][1] + rgb[1] * srgb2cieXYZ[2][2];
		// BP22 XYZ
		p22[2][0] = rgb[2] * srgb2cieXYZ[0][0] + rgb[2] * srgb2cieXYZ[0][1] + rgb[2] * srgb2cieXYZ[0][2];
		p22[2][1] = rgb[2] * srgb2cieXYZ[1][0] + rgb[2] * srgb2cieXYZ[1][1] + rgb[2] * srgb2cieXYZ[1][2];
		p22[2][2] = rgb[2] * srgb2cieXYZ[2][0] + rgb[2] * srgb2cieXYZ[2][1] + rgb[2] * srgb2cieXYZ[2][2];
	}

	std::vector<std::vector<double>> p22_out(3, std::vector<double>(3, 0.0));
	// transfer RGB gamut primaries to P22 gamut primaries while preserving luminosity
	if (direct_transfer) {
		p22_out = {
			{(RP22Chromaticity[0]),
			(RP22Chromaticity[1]),
			(RP22Chromaticity[2] * p22[0][1])},

			{(GP22Chromaticity[0]),
			(GP22Chromaticity[1]),
			(GP22Chromaticity[2] * p22[1][1])},

			{(BP22Chromaticity[0]),
			(BP22Chromaticity[1]),
			(BP22Chromaticity[2] * p22[2][1])}
		};
	}
	else {
		// convert to xyY
		// https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_xy_chromaticity_diagram_and_the_CIE_xyY_color_space
		p22_out = {
			// RP22 xyY
			{(p22[0][0] / (p22[0][0] + p22[0][1] + p22[0][2])),
			(p22[0][1] / (p22[0][1] + p22[0][1] + p22[0][2])),
			p22[0][1]},

			// GP22 xyY
			{(p22[1][0] / (p22[1][0] + p22[1][1] + p22[1][2])),
			(p22[1][1] / (p22[1][1] + p22[1][1] + p22[1][2])),
			p22[1][1]},

			// BP22 xyY
			{(p22[2][0] / (p22[2][0] + p22[2][1] + p22[2][2])),
			(p22[2][1] / (p22[2][1] + p22[2][1] + p22[2][2])),
			p22[2][1]}
		};
	}

	return p22_out;
}

std::vector<uint16_t> PhosphorFrame::P22TosRGB(std::vector<double> pixel_r, std::vector<double> pixel_g, std::vector<double> pixel_b, bool direct_transfer)
{
	// CIE xyY to CIE XYZ to sRGB
	// https://en.wikipedia.org/wiki/SRGB#From_CIE_XYZ_to_sRGB
	
	// convert to XYZ
	// https://en.wikipedia.org/wiki/CIE_1931_color_space#CIE_xy_chromaticity_diagram_and_the_CIE_xyY_color_space
	double pixel_XYZ[3][3] = {
		{(pixel_r[2] / pixel_r[1]) * pixel_r[0],
		pixel_r[2],
		(pixel_r[2] / pixel_r[1]) * (1 - pixel_r[0] - pixel_r[1])},

		{(pixel_g[2] / pixel_g[1]) * pixel_g[0],
		pixel_g[2],
		(pixel_g[2] / pixel_g[1]) * (1 - pixel_g[0] - pixel_g[1])},

		{(pixel_b[2] / pixel_b[1]) * pixel_b[0],
		pixel_b[2],
		(pixel_b[2] / pixel_b[1]) * (1 - pixel_b[0] - pixel_b[1])},
	};

	// clamp XYZ
	for (int i = 0; i < 3; i++) {
		if (pixel_XYZ[i][0] >= 0.9505) pixel_XYZ[i][0] = 0.9505;
		if (pixel_XYZ[i][1] >= 1.0) pixel_XYZ[i][1] = 1.0;
		if (pixel_XYZ[i][2] >= 1.0890) pixel_XYZ[i][2] = 1.0890;
	}

	double rgb[3]{};
	// convert P22 gamut to sRGB gamut
	if (direct_transfer) {
		// convert to RGB based on luminosity of primaries
		rgb[0] += pixel_r[2];
		rgb[1] += pixel_g[2];
		rgb[2] += pixel_b[2];
	}
	else {
		// CIE XYZ to linear RGB via matrix multiplication
		double ciexyz2srgb[3][3] = {
			{3.2406, -1.5372, -0.4986},
			{-0.9689, 1.8758, 0.0415},
			{0.05557, -0.2040, 1.0570}
		};

		rgb[0] += (pixel_XYZ[0][0] * ciexyz2srgb[0][0] + pixel_XYZ[0][1] * ciexyz2srgb[0][1] + pixel_XYZ[0][2] * ciexyz2srgb[0][2]);
		rgb[0] += (pixel_XYZ[1][0] * ciexyz2srgb[0][0] + pixel_XYZ[1][1] * ciexyz2srgb[0][1] + pixel_XYZ[1][2] * ciexyz2srgb[0][2]);
		rgb[0] += (pixel_XYZ[2][0] * ciexyz2srgb[0][0] + pixel_XYZ[2][1] * ciexyz2srgb[0][1] + pixel_XYZ[2][2] * ciexyz2srgb[0][2]);

		rgb[1] += (pixel_XYZ[0][0] * ciexyz2srgb[1][0] + pixel_XYZ[0][1] * ciexyz2srgb[1][1] + pixel_XYZ[0][2] * ciexyz2srgb[1][2]);
		rgb[1] += (pixel_XYZ[1][0] * ciexyz2srgb[1][0] + pixel_XYZ[1][1] * ciexyz2srgb[1][1] + pixel_XYZ[1][2] * ciexyz2srgb[1][2]);
		rgb[1] += (pixel_XYZ[2][0] * ciexyz2srgb[1][0] + pixel_XYZ[2][1] * ciexyz2srgb[1][1] + pixel_XYZ[2][2] * ciexyz2srgb[1][2]);

		rgb[2] += (pixel_XYZ[0][0] * ciexyz2srgb[2][0] + pixel_XYZ[0][1] * ciexyz2srgb[2][1] + pixel_XYZ[0][2] * ciexyz2srgb[2][2]);
		rgb[2] += (pixel_XYZ[1][0] * ciexyz2srgb[2][0] + pixel_XYZ[1][1] * ciexyz2srgb[2][1] + pixel_XYZ[1][2] * ciexyz2srgb[2][2]);
		rgb[2] += (pixel_XYZ[2][0] * ciexyz2srgb[2][0] + pixel_XYZ[2][1] * ciexyz2srgb[2][1] + pixel_XYZ[2][2] * ciexyz2srgb[2][2]);
	}

	// clamp RGB
	rgb[0] = rgb[0] <= 1.0 ? rgb[0] : 1.0;
	rgb[0] = rgb[0] >= 0.0 ? rgb[0] : 0.0;
	rgb[1] = rgb[1] <= 1.0 ? rgb[1] : 1.0;
	rgb[1] = rgb[1] >= 0.0 ? rgb[1] : 0.0;
	rgb[2] = rgb[2] <= 1.0 ? rgb[2] : 1.0;
	rgb[2] = rgb[2] >= 0.0 ? rgb[2] : 0.0;

	// linear RGB to sRGB
	for (int i = 0; i < 3; i++) {
		if (rgb[i] <= 0.0031308) rgb[i] *= 12.92;
		else rgb[i] = (std::pow(rgb[i], (1.0 / 2.4))) * 1.055 - 0.055;
	}

	uint16_t out_r = static_cast<uint16_t>(rgb[0] * static_cast<double>(0xFFFF));
	uint16_t out_g = static_cast<uint16_t>(rgb[1] * static_cast<double>(0xFFFF));
	uint16_t out_b = static_cast<uint16_t>(rgb[2] * static_cast<double>(0xFFFF));

	return { out_r, out_g, out_b };
}
