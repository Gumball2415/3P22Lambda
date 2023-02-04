#pragma once

class PhosphorFrame {
private:
	uint32_t Width = 0;
	uint32_t Height = 0;
	uint32_t Bitdepth = 0;
	double RP22Falloff = 0;
	double GP22Falloff = 0;
	double BP22Falloff = 0;
	std::vector<double> RP22FrameBuffer;
	std::vector<double> GP22FrameBuffer;
	std::vector<double> BP22FrameBuffer;
	std::vector<double> RP22FrameBufferPrev;
	std::vector<double> GP22FrameBufferPrev;
	std::vector<double> BP22FrameBufferPrev;
	std::vector<uint8_t> FrameBufferFinal;
	double ProcessPixel(double input_pixel, double phosphor_falloff, double previous_pixel);
public:
	PhosphorFrame(uint32_t frame_width, uint32_t frame_height, uint32_t frame_bitdepth, double rp22_falloff, double gp22_falloff, double bp22_falloff);
	~PhosphorFrame();
	void ProcessFrame(std::string file_path, std::string output_directory);
};
