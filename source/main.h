#pragma once

class PhosphorFrame {
private:
	uint32_t Width = 0;
	uint32_t Height = 0;
	uint32_t Bitdepth = 0;
	double RP22Falloff = 0.0;
	double GP22Falloff = 0.0;
	double BP22Falloff = 0.0;
	double FPSDelta = 0.0;

	// monitor primaries, in CIE xyY
	// based on RP 145-2004
	double RP22Chromaticity[3] = { 0.630, 0.340, 1.0 };
	double GP22Chromaticity[3] = { 0.310, 0.595, 1.0 };
	double BP22Chromaticity[3] = { 0.155, 0.070, 1.0 };

	std::vector<std::vector<double>> RP22FrameBuffer;
	std::vector<std::vector<double>> GP22FrameBuffer;
	std::vector<std::vector<double>> BP22FrameBuffer;
	std::vector<uint8_t> FrameBufferFinal;

	// input / output pixel range is 0-1
	std::vector<double> ProcessPixel(std::vector<double> input_pixel, double phosphor_falloff, std::vector<double> previous_pixel);

	// converts sRGB pixels to equivalent P22 color in CIE xyY coordinates
	std::vector<std::vector<double>> sRGBToP22(uint16_t pixel_r, uint16_t pixel_g, uint16_t pixel_b, bool direct_transfer);

	// converts P22 pixels in CIE xyY coordinates to equivalent color in sRGB colorspace
	std::vector<uint16_t> P22TosRGB(std::vector<double> pixel_r, std::vector<double> pixel_g, std::vector<double> pixel_b, bool direct_transfer);
public:
	PhosphorFrame(uint32_t frame_width, uint32_t frame_height, double rp22_falloff, double gp22_falloff, double bp22_falloff, double fps_delta);
	~PhosphorFrame();

	void ProcessFrame(std::string file_path, std::string output_directory, bool srgb_direct_transfer, bool p22_direct_transfer);
};
