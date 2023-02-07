#pragma once

class PhosphorFrame {
private:
	uint32_t Width = 0;
	uint32_t Height = 0;
	uint32_t Bitdepth = 0;

	// fade time, in milliseconds
	double RP22Falloff = 0.0;
	double GP22Falloff = 0.0;
	double BP22Falloff = 0.0;

	// FPS to calculate fade delta
	double FPS = 0.0;

	std::vector<std::vector<double>> P22ShortFrameBuffer;
	std::vector<std::vector<double>> P22LongFrameBuffer;
	std::vector<uint8_t> FrameBufferFinal;

	// input / output pixel range is 0-1
	std::vector<double> ProcessPixel(std::vector<double> in, std::vector<double> prev, double rp22_falloff, double gp22_falloff, double bp22_falloff);

	// converts sRGB pixels to three floating point linear RGB
	std::vector<double> sRGBToP22(uint16_t pixel_r, uint16_t pixel_g, uint16_t pixel_b);

	// converts P22 pixels in CIE xyY coordinates to equivalent color in sRGB colorspace
	std::vector<uint16_t> P22TosRGB(double pixel_r, double pixel_g, double pixel_b);
public:
	PhosphorFrame(uint32_t frame_width, uint32_t frame_height, double fps_delta);
	~PhosphorFrame();

	void ProcessFrame(std::string file_path, std::string output_directory);
};
