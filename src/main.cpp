#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <vector>
#include <math.h>
#include <iostream>
#include <string>

const double sigma_s = 9.0;		// \in [3, 15]
const double sigma_r = 0.07;	// \in [0.02, 0.09]
const double alpha = 0.5;

struct Image {
	int width, height, channels;
	std::vector<std::vector<std::vector<double>>> data; // [0.0, 1.0] for each value

	/***
	* constructor
	*/
	Image() {}

	/***
	* constructor
	* @param _width image width
	* @param _height image height
	* @param _channels image channel 
	*/
	Image(int _width, int _height, int _channels) {
		width = _width;
		height = _height;
		channels = _channels;
		data = std::vector<std::vector<std::vector<double>>>(
			height, std::vector<std::vector<double>>(
				width, std::vector<double>(
					channels, 0.0)));
	}

	/***
	* read .png image file
	* @param filename file name of the .png image
	*/
	void read_image(const char *filename) {
		unsigned char *image_data = stbi_load(filename, &width, &height, &channels, 0);
		if (!image_data) std::cerr << "Failed to load image data" << std::endl;

		data = std::vector<std::vector<std::vector<double>>>(
			height, std::vector<std::vector<double>>(
				width, std::vector<double>(
					channels, 0.0)));
		for (int h = 0; h < height; h++) {
			for (int w = 0; w < width; w++) {
				int index = (h*width + w) * channels;
				for (int c = 0; c < channels; c++) 
					data[h][w][c] = static_cast<double>(image_data[index + c]) / 255.0;
			}
		}
		stbi_image_free(image_data);
	}

	/***
	* write the current state to .png image file
	* @param filename filename of the .png image
	*/
	void write_image(const char *filename) {
		std::vector<unsigned char> image_data(height*width*channels);
		for (int h = 0; h < height; h++) {
			for (int w = 0; w < width; w++) {
				int index = (h*width + w) * channels;
				for (int c = 0; c < channels; c++) {
					unsigned char pixel_value = static_cast<unsigned char>(data[h][w][c] * 255.0);
					image_data[index + c] = std::max((unsigned char)(0), std::min((unsigned char)255, pixel_value));
				}
			}
		}

		int ok = stbi_write_png(filename, width, height, channels, image_data.data(), width*channels);
		if (!ok) std::cerr << "Failed to write image data" << std::endl;
	}
};

/*** 
* Image addition
* @param I left hand
* @param J right hand
* @return I + J
*/
Image image_add(Image &I, Image &J) {
	assert(I.width == J.width);
	assert(I.height == J.height);

	Image I_out(I.width, I.height, I.channels);
	for (int h = 0; h < I.width; h++) {
		for (int w = 0; w < I.height; w++) {
			for (int c = 0; c < I.channels; c++) {
				I_out.data[h][w][c] = I.data[h][w][c] + J.data[h][w][c];
				I_out.data[h][w][c] = std::min(1.0, std::max(0.0, I_out.data[h][w][c]));
			}
		}
	}
	return I_out;
}

/*** 
* Image subtraction
* @param I left hand
* @param J right hand
* @return I - J
*/
Image iamge_sub(Image &I, Image &J) {
	assert(I.width == J.width);
	assert(I.height == J.height);

	Image I_out(I.width, I.height, I.channels);
	for (int h = 0; h < I.width; h++) {
		for (int w = 0; w < I.height; w++) {
			int index = (h*I.width + w) * I.channels;
			for (int c = 0; c < I.channels; c++) {
				I_out.data[h][w][c] = I.data[h][w][c] - J.data[h][w][c];
				I_out.data[h][w][c] = std::min(1.0, std::max(0.0, I_out.data[h][w][c]));
			}
		}
	}
	return I_out;
}

/*** 
* return gaussian kernel in pixel space
* @param kernel_size gaussian kernel size
* @param sigma standard deviation
* @return gaussian kernel sized kernel_size
*/
std::vector<std::vector<double>> compute_gaussian_kernel(
	int kernel_size,
	double sigma
) {
	assert(kernel_size % 2 == 1);
	std::vector<std::vector<double>> kernel(kernel_size, std::vector<double>(kernel_size, 0.0));
	for (int h = 0; h < kernel_size; h++) {
		for (int w = 0; w < kernel_size; w++) {
			double dx = h - kernel_size/2.0;
			double dy = w - kernel_size/2.0;
			double dist = dx*dx + dy*dy;

			kernel[h][w] = std::exp(-dist / (2.0*sigma*sigma));
		}
	}

	return kernel;
}

/***
* apply bilinear interpolation
* @param I_in input image
* @param I_out output imgae
* @param alpha magnification power
*/
void bilinear_interpolation(
	Image &I_in,
	Image &I_out,
	double alpha
) {
	assert(static_cast<int>(I_in.width*alpha) == I_out.width);
	assert(static_cast<int>(I_in.height*alpha) == I_out.height);

	for (int h = 0; h < I_out.height; h++) {
		for (int w = 0; w < I_out.width; w++) {
			for (int c = 0; c < I_out.channels; c++) {
				int h_ = floor(h / alpha);
				int w_ = floor(w / alpha);
				double dh = h / alpha - static_cast<double>(h_);
				double dw = w / alpha - static_cast<double>(w_);

				if (h_ + 1 < I_in.height && w_ + 1 < I_in.width) {
					I_out.data[h][w][c] += (1.0 - dh)*(1.0 - dw)*I_in.data[h_][w_][c];
					I_out.data[h][w][c] += dh*(1.0 - dw)*I_in.data[h_ + 1][w_][c];
					I_out.data[h][w][c] += (1.0 - dh)*dw*I_in.data[h_][w_ + 1][c];
					I_out.data[h][w][c] += dh*dw*I_in.data[h_ + 1][w_ + 1][c];
				} else if (h_ + 1 < I_in.height) {
					I_out.data[h][w][c] += (1.0 - dh)*I_in.data[h_][w_][c];
					I_out.data[h][w][c] += dh*I_in.data[h_ + 1][w_][c];
				} else if (w_ + 1 < I_in.width) {
					I_out.data[h][w][c] += (1.0 - dw)*I_in.data[h_][w_][c];
					I_out.data[h][w][c] += dw*I_in.data[h_][w_ + 1][c];
				} else {
					I_out.data[h][w][c] = I_in.data[h_][w_][c];
				}

				I_out.data[h][w][c] = std::min(1.0, std::max(0.0, I_out.data[h][w][c]));
			}
		}
	}
}

/***
* apply gaussian smoothing
* @param I_in input image
* @param I_out output image
* @param kernel_size gaussian kernel size
* @param sigma standard deviation
*/
void gaussian_smoothing(
	Image &I_in,
	Image &I_out, 
	int kernel_size, 
	double sigma
) {
	assert(I_in.width == I_out.width);
	assert(I_in.height == I_out.height);

	// precomputed gaussian kernel
	std::vector<std::vector<double>> gaussian_kernel = compute_gaussian_kernel(kernel_size, sigma);

	// gaussian smoothing
	for (int h = 0; h < I_in.height; h++) {
		for (int w = 0; w < I_in.width; w++) {
			for (int c = 0; c < I_in.channels; c++) {
				double g_sum = 0.0;
				for (int i = -kernel_size/2; i < kernel_size/2; i++) {
					for (int j = -kernel_size/2; j < kernel_size/2; j++) {
						int h_ = h + i;
						int w_ = w + j;
						if (h_ < 0 || h_ >= I_in.height || w_ < 0 || w_ >= I_in.width) continue;

						double g = gaussian_kernel[i + kernel_size/2][j + kernel_size/2];
						I_out.data[h][w][c] += g * I_in.data[h_][w_][c];
						g_sum += g;
					}
				}

				// normalization
				assert(g_sum != 0.0);
				I_out.data[h][w][c] /= g_sum;
				I_out.data[h][w][c] = std::min(1.0, std::max(0.0, I_out.data[h][w][c]));
			}
		}
	}
}

/*** 
* apply joint bilateral filter. switch between joint bilateral filter and 
* joint bilateral upsampling using alpha.
* @param I input image
* @param J guidance image
* @param S output image
* @param kernel_size gaussian kernel size
* @param sigma_s standard deviation for position distance
* @param sigma_r standard deviation for pixel distance
*/
void joint_bilateral_filter(
	Image &I, 
	Image &J, 
	Image &S,
	int kernel_size,
	double sigma_s,
	double sigma_r
) {
	assert(I.width == J.width);
	assert(I.height == J.height);
	assert(J.width == S.width);
	assert(J.height == S.height);

	// precomputed gaussian kernel
	std::vector<std::vector<double>> gaussian_kernel = compute_gaussian_kernel(kernel_size, sigma_s);

	for (int h = 0; h < I.height; h++) {
		for (int w = 0; w < I.width; w++) {
			double K_p = 0.0;
			for (int i = -kernel_size/2; i < kernel_size/2; i++) {
				for (int j = -kernel_size/2; j < kernel_size/2; j++) {
					int h_ = h + i;
					int w_ = w + j;
					if (h_ < 0 || h_ >= I.height || w_ < 0 || w_ >= I.width) continue;

					// compute g_s
					double g_s = gaussian_kernel[i + kernel_size/2][j + kernel_size/2];

					// compute g_r
					double I_dist = 0.0;
					for (int c = 0; c < I.channels; c++) {
						double dc = J.data[h][w][c] - J.data[h_][w_][c];
						I_dist += dc * dc;
					}
					double g_r = std::exp(-I_dist / (2.0*sigma_r*sigma_r));

					// add the computed value to S[h][w]
					for (int c = 0; c < I.channels; c++) S.data[h][w][c] += g_s*g_r*I.data[h_][w_][c];

					// add g_s*g_r to K_p
					K_p += g_s*g_r;
				}
			}
			
			// normalization
			assert(K_p != 0.0);
			for (int c = 0; c < I.channels; c++) {
				S.data[h][w][c] /= K_p;
				S.data[h][w][c] = std::min(1.0, std::max(0.0, S.data[h][w][c]));
			}
		}
	}
}

/*** 
* downsample the input image
* @param I_in input image
* @param I_out downsampled image
* @param alpha magnification power
*/
void downsample(Image &I_in, Image &I_out, double alpha) {
	// 1. gaussian smoothing
	int kernel_size = 5;
	double sigma = 1.0;
	Image I_gaus(I_in.width, I_in.height, I_in.channels);
	gaussian_smoothing(I_in, I_gaus, kernel_size, sigma);

	// 2. downsampling with bilinear interpolation
	bilinear_interpolation(I_gaus, I_out, alpha);
}

/*** 
* upsample the input image
* @param I_in input image
* @param I_out upsampled image
* @param alpha magnification power
*/
void upsample(Image &I_in, Image &I_out, double alpha) {
	// 1. upsampling with bilinear interpolation
	Image I_up(I_in.width*alpha, I_in.height*alpha, I_in.channels);
	bilinear_interpolation(I_in, I_up, alpha);

	// 2. gaussian smoothing
	int kernel_size = 5;
	double sigma = 1.0;
	gaussian_smoothing(I_up, I_out, kernel_size, sigma);
}

/*** 
* debug function for bilinear_interpolation()
* @param G_0 input image
*/
void bilinear_interpolation_debug(Image &G_0, double alpha) {
	double down_alpha = alpha;
	Image G_down(G_0.height*down_alpha, G_0.width*down_alpha, G_0.channels);
	bilinear_interpolation(G_0, G_down, down_alpha);
	const char *G_down_filename = "./debug/G_down.png";
	G_down.write_image(G_down_filename);

	double up_alpha = 1.0 / alpha;
	Image G_up(G_down.height*up_alpha, G_down.width*up_alpha, G_0.channels);
	bilinear_interpolation(G_down, G_up, up_alpha);
	const char *G_up_filename = "./debug/G_up.png";
	G_up.write_image(G_up_filename);

	Image L = iamge_sub(G_0, G_up);
	const char *L_filename = "./debug/L.png";
	L.write_image(L_filename);
}

/*** 
* debug function for gaussian_smoothing()
* @param G_0 input image
*/
void gaussian_smoothing_debug(Image &G_0) {
	int kernel_size = 5;
	double sigma = 1.0;
	Image G_gaus(G_0.width, G_0.height, G_0.channels);
	gaussian_smoothing(G_0, G_gaus, kernel_size, sigma);

	const char *G_gaus_filename = "./debug/G_gaus.png";
	G_gaus.write_image(G_gaus_filename);
}

/***
* debug function for Gaussian/Laplacian pyramid
* @param G Gaussian pyramid
* @param L Laplacian pyramid
* @param N pyramid height
*/
void gaussian_laplacian_pyramid_debug(
	std::vector<Image> &G, 
	std::vector<Image> &L, 
	int N
) {	
	std::string G_filename = "./debug/G_";
	for (int k = 0; k < N + 1; k++) {
		std::string G_k_filename = G_filename + std::to_string(k) + ".png";
		G[k].write_image(G_k_filename.c_str());
	}

	std::string L_filename = "./debug/L_";
	for (int k = 0; k < N + 1; k++) {
		std::string L_k_filename = L_filename + std::to_string(k) + ".png";
		L[k].write_image(L_k_filename.c_str());
	}
}

/*** 
* debug function for joint_bilateral_filter()
* @param R_up upsampled image
* @param R_hat JBF(R_up + G)
* @param R_inter R_hat + L
* @param R JBF(R_hat + L, R_hat)
* @param N pyramid height
*/
void joing_bilateral_filter_debug(
	std::vector<Image> &R_up,
	std::vector<Image> &R_hat,
	std::vector<Image> &R_inter, 
	std::vector<Image> &R, 
	int N
) {
	std::string R_up_filename = "./debug/R_up_";
	for (int k = 0; k < N; k++) {
		std::string R_up_k_filename = R_up_filename + std::to_string(k) + ".png";
		R_up[k].write_image(R_up_k_filename.c_str());
	}

	std::string R_hat_filename = "./debug/R_hat_";
	for (int k = 0; k < N; k++) {
		std::string R_hat_k_filename = R_hat_filename + std::to_string(k) + ".png";
		R_hat[k].write_image(R_hat_k_filename.c_str());
	}

	std::string R_inter_filename = "./debug/R_inter_";
	for (int k = 0; k < N; k++) {
		std::string R_inter_k_filename = R_inter_filename + std::to_string(k) + ".png";
		R_inter[k].write_image(R_inter_k_filename.c_str());
	}
	
	std::string R_filename = "./debug/R_";
	for (int k = 0; k < N; k++) {
		std::string R_k_filename = R_filename + std::to_string(k) + ".png";
		R[k].write_image(R_k_filename.c_str());
	}
}

/*** 
* log the refined image
* @param R refined image
* @param filename original file name
*/
void log(Image &R, const char *filename) {
	std::string filename_str(filename);
	int last_delim = -1;
	for (int i = 0; i < filename_str.size(); i++) {
		if (filename_str[i] == '/') last_delim = i;
	}

	std::string R_filename = std::string("./output") + filename_str.substr(last_delim);
	R.write_image(R_filename.c_str());
	std::cout << "exported " << R_filename << std::endl;
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "Please pass the input file name." << std::endl;
		return 1;
	}
	const char *filename = argv[1];

	// Instantiate the input image
	Image G_0; 
	G_0.read_image(filename);

	/*** 
	* 0. determine the height of the gaussian & laplacian pyramid
	*		limit the long axis of the coarsest Gaussian pyramid level
	*		to [32, 64)
	*/
	int N = 0;
	int curr_width = G_0.width;
	int curr_height = G_0.height;
	while (curr_width >= 64 || curr_height >= 64) {
		curr_width /= 2;
		curr_height /= 2;
		N += 1;
	}

	/***
	*	1. construct Gaussian & Laplacian pyramid
	*/
	std::vector<Image> G(N + 1);
	G[0] = G_0;
	// initialization of G[1] ... G[N]
	for (int k = 0; k < N; k++) {
		std::cout << "Gaussian pyramid: " << k << std::endl;

		double down_alpha = 0.5;
		G[k + 1] = Image(G[k].width*down_alpha, G[k].height*down_alpha, G[k].channels);
		downsample(G[k], G[k + 1], down_alpha);
	}

	std::vector<Image> L(N + 1);
	L[N] = G[N];
	// initialization of L[0] ... L[N - 1]
	for (int k = 0; k < N; k++) {
		std::cout << "Laplacian pyramid: " << k << std::endl;

		double up_alpha = 2.0;
		Image upsampled_G = Image(G[k + 1].width*up_alpha, G[k + 1].width*up_alpha, G[k + 1].channels);
		upsample(G[k + 1], upsampled_G, up_alpha);
		L[k] = iamge_sub(G[k], upsampled_G);
	}

	/*** 
	* 2. upsampling for texture filitering
	*/
	std::vector<Image> R_up(N + 1);
	std::vector<Image> R_hat(N + 1);
	std::vector<Image> R_inter(N + 1);
	std::vector<Image> R(N + 1);

	R[N] = G[N];
	for (int k = N; k > 0; k--) {
		std::cout << "PSU: " << k << std::endl;

		double sigma_s_k = sigma_s / std::pow(2.0, k - 1);
		int kernel_size_1 = static_cast<int>(std::round(std::max(sigma_s_k, 3.0)));
		int kernel_size_2 = static_cast<int>(std::round(std::max(4.0*sigma_s_k, 3.0)));
		if (kernel_size_1 % 2 == 0) {
			if (std::max(sigma_s_k, 3.0) > kernel_size_1) kernel_size_1 += 1;
			else kernel_size_1 -= 1;
		}
		if (kernel_size_2 % 2 == 0) {
			if (std::max(4.0*sigma_s_k, 3.0) > kernel_size_2) kernel_size_2 += 1;
			else kernel_size_2 -= 1;
		}
		double up_alpha = 2.0;

		// joint bilateral upsampling
		R_up[k - 1] = Image(R[k].width*up_alpha, R[k].height*up_alpha, R[k].channels);
		R_hat[k - 1] = Image(R[k].width*up_alpha, R[k].height*up_alpha, R[k].channels);
		bilinear_interpolation(R[k], R_up[k - 1], up_alpha);
		joint_bilateral_filter(R_up[k - 1], G[k - 1], R_hat[k - 1], kernel_size_1, sigma_s_k, sigma_r);

		// joint bilateral filitering
		R[k - 1] = Image(R[k].width*up_alpha, R[k].height*up_alpha, R[k].channels);
		R_inter[k - 1] = image_add(R_hat[k - 1], L[k - 1]);
		joint_bilateral_filter(R_inter[k - 1], R_hat[k - 1], R[k - 1], kernel_size_2, sigma_s_k, sigma_r);

		// enhancement
		Image R_ref = Image(R[k].width*up_alpha, R[k].height*up_alpha, R[k].channels);
		joint_bilateral_filter(R[k - 1], R[k - 1], R_ref, kernel_size_2, sigma_s_k, sigma_r);
		R[k - 1] = R_ref;
	}

	/*** 
	* 3. debug & log
	*/
	log(R[0], filename);
	bilinear_interpolation_debug(G_0, alpha);
	gaussian_smoothing_debug(G_0);
	gaussian_laplacian_pyramid_debug(G, L, N);
	joing_bilateral_filter_debug(R_up, R_hat, R_inter, R, N);
}