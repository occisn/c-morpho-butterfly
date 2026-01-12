#include <assert.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define IMAGE_HEIGHT 1100 // for n
#define IMAGE_WIDTH 2000  // for m
#define OUTPUT_DIR "pics"

static inline double exp_minus_exp(double x)
{
  return (x > 85.0) ? 0.0 : exp(-exp(x));
}

static inline double exp_minus_exp_minus_exp(double x, double y)
{
  return (x > 85.0 || y > 85.0) ? 0.0 : exp(-exp(x) - exp(y));
}

/* Calculates C(x,y) */
static inline double C(double x, double y)
{
  const double res = pow(sin(0 + (14 * atan((100 * (y + (x * 0.25) - (1.0 / 25.0))) / (1.0 + fabs(100.0 * x - 25.0 * y - 3.0 * atan(100.0 * x - 25.0 * y))))) + (14 * fabs((x * 0.5) - (y / 8.0)))), 4);
  return res;
}

/* Calculates E(x,y) */
static inline double E(double x, double y)
{
  const double res = 1.0 - exp(-exp(0 + (-100.0 * pow((3.0 * y + 0.75 * x + 0.27), 4)) + (-100.0 * pow(fabs(7 * (1.0 + 1.0 / (sqrt(fabs(100.0 * y + 25.0 * x - 6.0)) + 0.3)) * (x - y * 0.25)), (3.0 * y + 0.75 * x + 2.27))) + 10.0)) * (1.0 - exp_minus_exp_minus_exp((200.0 * fabs(y + x * 0.25 - 0.2 + 3 * (x - y * 0.25) * (x - y * 0.25)) - 32.0), (500.0 * fabs(y + x * 0.25 - (1.0 / 20.0) - (0.7 * sqrt(fabs(x - y * 0.25)))) - 2.5)));
  return res;
}

/* Calculates L(x,y) */
static inline double L(double x, double y)
{
  double res = 0.0;
  for (uint8_t s = 1; s <= 25; s++) {
    res += pow(sin(((80.0 + 30.0 * sin(s * s)) * atan((100.0 * y + 25.0 * x - 4.0 * sin(s)) / (1.0 + fabs(100.0 * x - 25.0 * y - 3.0 * atan(100.0 * x - 25.0 * y))))) + fabs((x * 0.5) - (y / 8.0)) + (4.0 * sin(5.0 * s))), 6);
  }
  return res;
}

/* Calculates W(x,y) knowing C(x,y) */
static inline double W(double x, double y, double Cxy)
{
  const double omega1 = 0.0 + (-40 * Cxy) + (196.0 / 5.0) + ((4.0 / 5.0) * (sqrt(0 + ((x - (y * 0.25)) * (x - (y * 0.25))) + ((y + (x * 0.25)) * (y + (x * 0.25))))));
  const double omega2 = -(40.0 * (0 + (5.0 * fabs(y + (x * 0.25) - (3.0 / 50.0) + ((1.0 / 3.0) * (x - (y * 0.25)) * (x - (y * 0.25))))) + pow(fabs((2.0 * x) - (y * 0.5)), 3) - (2.0 / 5.0)));
  const double omega3 = -1000 * fabs(x - (y * 0.25)) + 100.0 - 90.0 * atan(8.0 * y + 2.0 * x + (8.0 / 5.0));
  const double omega4 = 1000 * (0.0 + fabs(x - y * 0.25) - (7.0 / 50.0) + (9.0 * (y + (x * 0.25) + 0.2) * (1.0 / 20.0)));
  const double omega5 = 70 * fabs((5 * (fabs(y + x * 0.25 - (3.0 / 50.0) + ((1.0 / 3.0) * (x - (y * 0.25)) * (x - (y * 0.25)))))) + (pow(fabs((2.0 * x) - (y * 0.5)), 3)) - (2.0 / 5.0)) - (1.0 / 200.0);
  const double omega6 = 700 * fabs(0.0 + fabs(x - y * 0.25) - 0.1 + (0.9 * atan(8.0 * (y + (x * 0.25) + 0.2)))) - (21.0 / 20.0);
  const double res = (-exp_minus_exp_minus_exp(omega1, omega2)) * (1.0 - exp_minus_exp_minus_exp(omega3, omega4)) - exp_minus_exp(omega5) - exp_minus_exp(omega6) + 1.0;
  return res;
}

/* Calculates A(v,x,y) knowing C(x,y)*/
static inline void A(double *Axy, double x, double y, double Cxy)
{
  const double A_part1 = (x * 0.25) + y - 0.25 * fabs(sin((12.0 / 5.0) * ((0.7 * fabs(x - (y * 0.25)) + (0.3 * sqrt(fabs(x - y * 0.25)))))));
  const double A_part2 = x * 0.25 + (7.0 / 20.0) + y + 0.2 * atan(6.0 * fabs(x - (y * 0.25))) + 0.2 * atan(40.0 * fabs(x - (y * 0.25))) - (23.0 / 20.0) * (1.5 + (1.0 / 25.0) * cos(10.0 * (y + (x * 0.25) + (6.0 / 25.0))) + 0.03 * Cxy + 0.3 * atan(30.0 * (y + (x * 0.25) - 0.25))) * fabs(x - (y * 0.25));
  for (uint8_t v = 0; v <= 1; v++) {
    Axy[v] = exp(-exp(200.0 * (v * (1.0 / 50.0) + A_part1)) - exp(-200.0 * (A_part2 - v * (7.0 / 50.0))));
  }
}

/* Calculates K(v,x,y) */
static inline void K(double *Kxy, double x, double y)
{
  for (uint8_t v = 0; v <= 2; v++) {
    double Kvxy = 0.0;
    for (uint8_t i = 1; i <= 60; i++) {
      double s = i;
      Kvxy += ((5.0 / 2.0) * ((2.0 / 25.0) + (3.0 / 50.0) * cos(s * (4.0 + 4.0 * v))) * (sin(5.0 * s) + sin(2.0 * s) + 3.0) * 0.2 * exp(-exp(-25.0 * (pow(sin(sin(2.0 * s) + (6 + sin(s * s)) * (sin(7.0 * s) * (x * 0.5) + cos(7.0 * s) * ((y - 8.0) * 0.5))), 10) * pow(sin(sin(3.0 * s) + (6 + 2 * sin(s * s)) * (sin(7.0 * s) * ((y - 8) * 0.5) - cos(7.0 * s) * (x * 0.5))), 10) - 0.1))));
    }
    Kxy[v] = Kvxy;
  }
}

/* Calculates H(v,x,y) knowing E(x,y), L(x,y), W(x,y), A(v,x,y), K(v,x,y) */
static inline void H(double *Hxy, double x, double y, double Exy, double Lxy, double Wxy, const double *Axy, const double *Kxy)
{
  // const double H_part1 = Axy[0] * Axy[1] * (1.0 - Exy) * (1.0 + Lxy * (1.0 / 50.0)) * exp_minus_exp(0 + exp((2.0 * y) + (0.5 * x) + (2.0 / 5.0) - (2.0 * fabs(x - (y * 0.25)))) + exp((8.0 * y) + (2.0 * x) + (2.0 / 5.0) - fabs((8.0 * x) - (2.0 * y)))) * Wxy;
  const double H_part1 = Axy[0] * Axy[1] * (1.0 - Exy) * (1.0 + Lxy * (1.0 / 50.0)) * exp_minus_exp_minus_exp((2.0 * y) + (0.5 * x) + (2.0 / 5.0) - (2.0 * fabs(x - (y * 0.25))), (8.0 * y) + (2.0 * x) + (2.0 / 5.0) - fabs((8.0 * x) - (2.0 * y))) * Wxy;
  const double H_part2 = exp_minus_exp(-(50.0 *
                                         ((pow(cos((2.0 * y) + (x * 0.5) + (7.0 / 5.0) - fabs((2.0 * x) - (y * 0.5))),
                                               80) *
                                           pow(sin((20.0 * y) + (5.0 * x) + fabs((20.0 * x) - (5.0 * y))),
                                               2)) -
                                          pow(0.0 + (2.7 * y) + ((27.0 * x) * (1.0 / 40.0)) + (81.0 / 250.0),
                                              10) -
                                          (49.0 / 50.0))));
  for (uint8_t v = 0; v <= 2; v++) {
    double Kvxy = Kxy[v];
    Hxy[v] = 0.0 + (((18 - (9.0 * v) + (v * v)) * (1.0 / 20.0)) * (1.0 - Axy[0]) * (1.0 - Exy) * Kvxy) + ((2 + (3.0 * v)) * 0.2) * H_part1 + H_part2 + 0.1 * Exy * ((v - 1.0) * (v - 1));
  }
}

/* Calcules F(z) */
static inline double F(double z)
{
  const double res = floor(255.0 * exp_minus_exp(-(1000.0 * z)) * pow(fabs(z), exp_minus_exp(1000.0 * (z - 1.0))));
  return res;
}

/**
 * Helper to convert Hue to RGB component
 */
static double hue_to_rgb(double p, double q, double t)
{
  if (t < 0.0)
    t += 1.0;
  if (t > 1.0)
    t -= 1.0;
  if (t < 1.0 / 6.0)
    return p + (q - p) * 6.0 * t;
  if (t < 0.5)
    return q;
  if (t < 2.0 / 3.0)
    return p + (q - p) * (2.0 / 3.0 - t) * 6.0;
  return p;
}

/**
 * Heatmap Generator
 * @param height        Image height
 * @param width         Image width
 * @param value_array   A pointer to a array of height x width
 * @param export_file   Filename string
 */
bool draw_heatmap_from_values(int height, int width, double *value_array, const char *export_file)
{
  if (height <= 0 || width <= 0 || !value_array)
    return false;

  // (1) Calculate min, max, and range
  // We can now access elements using 2D syntax: value_array[row][col]
  double value_min = value_array[0];
  double value_max = value_array[0];

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      double z = value_array[y * width + x];
      if (z < value_min)
        value_min = z;
      if (z > value_max)
        value_max = z;
    }
  }

  double value_range = value_max - value_min;
  if (value_range == 0.0)
    value_range = 1.0;

  // (2) Calculate RGB values
  uint8_t *rgb_image = (uint8_t *)malloc(width * height * 3);
  if (!rgb_image)
    return false;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {

      double z = value_array[y * width + x];
      double w = (z - value_min) / value_range;

      double h = (w * 1.2) / 3.6;
      double s = 1.0;
      double l = 0.5;

      double q = (l < 0.5) ? (l * (1.0 + s)) : (l + s - l * s);
      double p = 2.0 * l - q;

      // Calculate flat index for the output buffer
      int pixel_idx = (y * width + x) * 3;
      rgb_image[pixel_idx + 0] = (uint8_t)(fmin(255.0, floor(255.0 * hue_to_rgb(p, q, h + 1.0 / 3.0))));
      rgb_image[pixel_idx + 1] = (uint8_t)(fmin(255.0, floor(255.0 * hue_to_rgb(p, q, h))));
      rgb_image[pixel_idx + 2] = (uint8_t)(fmin(255.0, floor(255.0 * hue_to_rgb(p, q, h - 1.0 / 3.0))));
    }
  }

  // (3) Export to JPG
  const int quality = 90;
  const int channels = 3;
  if (stbi_write_jpg(export_file, width, height, channels, rgb_image, quality) == 0) {
    perror("Problem in stbi_write_jpg\n");
    free(rgb_image);
    return false;
  }

  free(rgb_image);

  return true;
}

static bool generate_butterfly(double *duration)
{

  // (1) Start duration measurement

  struct timespec start = {0}, end = {0};
  clock_gettime(CLOCK_MONOTONIC, &start);

  // (2) Allocate RGB array

  uint8_t *rgb_array = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * 3);
  if (rgb_array == NULL) {
    perror("Problem in malloc of rgb_array\n");
    return false;
  }

  // (3) Populate RGB array

  const int num_cores = omp_get_num_procs();
  printf("\nNumber of cores on this system: %d\n", num_cores);

  printf("\nParallel populating of RGB arrays...\n");

  int chunk_size = IMAGE_HEIGHT / (4 * num_cores);
  if (chunk_size < 1) {
    chunk_size = 1;
  }

#pragma omp parallel for schedule(static, chunk_size)
  for (int n = 1; n <= IMAGE_HEIGHT; n++) {

    double Axy[2];
    double Hxy[3];
    double Kxy[3];

    for (int m = 1; m <= IMAGE_WIDTH; m++) {

      double x = (m - 1000.0) / 960.0;
      double y = (451.0 - n) / 960.0;
      double Cxy = C(x, y);
      double Exy = E(x, y);
      double Lxy = L(x, y);
      double Wxy = W(x, y, Cxy);
      A(Axy, x, y, Cxy);
      K(Kxy, x, y);
      H(Hxy, x, y, Exy, Lxy, Wxy, Axy, Kxy);

      int idx = ((n - 1) * IMAGE_WIDTH + (m - 1)) * 3;
      rgb_array[idx + 0] = F(Hxy[0]);
      rgb_array[idx + 1] = F(Hxy[1]);
      rgb_array[idx + 2] = F(Hxy[2]);

      if ((m == 1) && ((n == 1) || ((n % 100) == 0))) {
        printf("n = %d / %d\n", n, IMAGE_HEIGHT);
        fflush(stdout);
      }
    }
  }

  // (4) Create output image

  printf("Creating image...\n");

  const int channels = 3;
  const int quality = 90;
  if (stbi_write_jpg(OUTPUT_DIR "/output-butterfly.jpg",
                     IMAGE_WIDTH,
                     IMAGE_HEIGHT,
                     channels,
                     rgb_array,
                     quality) == 0) {
    perror("Problem in stbi_write_jpg\n");
    free(rgb_array);
    return false;
  }

  // (5) Free memory allocated to RGB array

  printf("Freeing memory allocated to RGB array...\n");
  free(rgb_array);

  // (6) Stop duration measurement

  clock_gettime(CLOCK_MONOTONIC, &end);
  *duration = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
  printf("Duration: %f seconds\n", *duration);
  fflush(stdout);

  return true;
}

static bool generate_heatmaps(void)
{

  // (1) Allocate heatmaps array

  double *C_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *C_values);
  if (C_values == NULL) {
    perror("Error in malloc of C_values\n");
    return false;
  }

  double *E_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *E_values);
  if (E_values == NULL) {
    perror("Error in malloc of E_values\n");
    return false;
  }

  double *L_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *L_values);
  if (L_values == NULL) {
    perror("Error in malloc of L_values\n");
    return false;
  }

  double *W_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *W_values);
  if (W_values == NULL) {
    perror("Error in malloc of W_values\n");
    return false;
  }

  double *A0_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *A0_values);
  if (A0_values == NULL) {
    perror("Error in malloc of A0_values\n");
    return false;
  }

  double *A1_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *A1_values);
  if (A1_values == NULL) {
    perror("Error in malloc of A1_values\n");
    return false;
  }

  double *K0_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *K0_values);
  if (K0_values == NULL) {
    perror("Error in malloc of K0_values\n");
    return false;
  }

  double *K1_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *K1_values);
  if (K1_values == NULL) {
    perror("Error in malloc of K1_values\n");
    return false;
  }

  double *K2_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *K2_values);
  if (K2_values == NULL) {
    perror("Error in malloc of K2_values\n");
    return false;
  }

  double *H0_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *H0_values);
  if (H0_values == NULL) {
    perror("Error in malloc of H0_values\n");
    return false;
  }

  double *H1_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *H1_values);
  if (H1_values == NULL) {
    perror("Error in malloc of H1_values\n");
    return false;
  }

  double *H2_values = malloc(IMAGE_HEIGHT * IMAGE_WIDTH * sizeof *H2_values);
  if (H2_values == NULL) {
    perror("Error in malloc of H2_values\n");
    return false;
  }

  // (2) Populate arrays

  const int num_cores = omp_get_num_procs();
  printf("\nNumber of cores on this system: %d\n", num_cores);

  printf("\nParallel calculations...\n");

  int chunk_size = IMAGE_HEIGHT / (4 * num_cores);
  if (chunk_size < 1) {
    chunk_size = 1;
  }

#pragma omp parallel for schedule(static, chunk_size)
  for (int n = 1; n <= IMAGE_HEIGHT; n++) {

    double Axy[2];
    double Kxy[3];
    double Hxy[3];

    for (int m = 1; m <= IMAGE_WIDTH; m++) {

      double x = (m - 1000.0) / 960.0;
      double y = (451.0 - n) / 960.0;
      double Cxy = C(x, y);
      C_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Cxy;
      double Exy = E(x, y);
      E_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Exy;
      double Lxy = L(x, y);
      L_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Lxy;
      double Wxy = W(x, y, Cxy);
      W_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Wxy;
      A(Axy, x, y, Cxy);
      A0_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Axy[0];
      A1_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Axy[1];
      K(Kxy, x, y);
      K0_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Kxy[0];
      K1_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Kxy[1];
      K2_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Kxy[2];
      H(Hxy, x, y, Exy, Lxy, Wxy, Axy, Kxy);
      H0_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Hxy[0];
      H1_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Hxy[1];
      H2_values[(n - 1) * IMAGE_WIDTH + (m - 1)] = Hxy[2];

      if ((m == 1) && ((n == 1) || ((n % 100) == 0))) {
        printf("n = %d / %d\n", n, IMAGE_HEIGHT);
        fflush(stdout);
      }
    }
  }

  // (3) Create heatmaps

  printf("\n");

  printf("Generate C heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, C_values, OUTPUT_DIR "/C-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for C heatmap\n");
    return false;
  }

  printf("Generate E heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, E_values, OUTPUT_DIR "/E-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for E heatmap\n");
    return false;
  }

  printf("Generate L heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, L_values, OUTPUT_DIR "/L-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for L heatmap\n");
    return false;
  }

  printf("Generate W heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, W_values, OUTPUT_DIR "/W-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for W heatmap\n");
    return false;
  }

  printf("Generate A0 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, A0_values, OUTPUT_DIR "/A0-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for A0 heatmap\n");
    return false;
  }

  printf("Generate A1 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, A1_values, OUTPUT_DIR "/A1-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for A1 heatmap\n");
    return false;
  }

  printf("Generate K0 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, K0_values, OUTPUT_DIR "/K0-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for K0 heatmap\n");
    return false;
  }

  printf("Generate K1 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, K1_values, OUTPUT_DIR "/K1-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for K1 heatmap\n");
    return false;
  }

  printf("Generate K2 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, K2_values, OUTPUT_DIR "/K2-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for K2 heatmap\n");
    return false;
  }

  printf("Generate H0 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, H0_values, OUTPUT_DIR "/H0-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for H0 heatmap\n");
    return false;
  }

  printf("Generate H1 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, H1_values, OUTPUT_DIR "/H1-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for H1 heatmap\n");
    return false;
  }

  printf("Generate H2 heatmap...\n");
  if (!draw_heatmap_from_values(IMAGE_HEIGHT, IMAGE_WIDTH, H2_values, OUTPUT_DIR "/H2-heatmap.jpg")) {
    perror("Problem in draw_heatmap_from_values for H2 heatmap\n");
    return false;
  }

  // (4) Free memory allocated to heatmaps

  printf("\nFreeing memory allocated to heatmaps...\n");
  free(C_values);
  free(E_values);
  free(L_values);
  free(W_values);
  free(A0_values);
  free(A1_values);
  free(K0_values);
  free(K1_values);
  free(K2_values);
  free(H0_values);
  free(H1_values);
  free(H2_values);

  return true;
}

bool run_butterfly(void)
{
  const int NB_RUNS = 1;
  double durations[NB_RUNS];
  double min_duration = DBL_MAX;

  for (int i = 0; i < NB_RUNS; i++) {
    printf("\nRun #%d\n", i + 1);
    printf("------\n");
    double duration0 = 0.0;
    if (!generate_butterfly(&duration0)) {
      perror("Problem in generate_butterfly");
      return false;
    }
    durations[i] = duration0;
    if (durations[i] < min_duration) {
      min_duration = durations[i];
    }
  }

  printf("\nAll execution durations:\n");
  for (int i = 0; i < NB_RUNS; i++) {
    printf("  Run %d: %f seconds\n", i + 1, durations[i]);
  }
  printf("Quickest execution duration: %f seconds\n", min_duration);
  return true;
}

int main(void)
{

  // (1) ensure output directory

  struct stat st = {0};
  if (stat(OUTPUT_DIR, &st) == -1) {
#ifdef _WIN32
    if (mkdir(OUTPUT_DIR) != 0) {
#else
    if (mkdir(OUTPUT_DIR, 0755) != 0) {
#endif
      perror("Failed to create output directory");
      return EXIT_FAILURE;
    }
    printf("Created output directory: %s\n", OUTPUT_DIR);
  }

  // (2) Generate butterfly

  printf("\n(1) GENERATE BUTTERFLY\n");
  printf("----------------------\n");

  if (!run_butterfly()) {
    perror("Problem in generating butterfly\n");
    return EXIT_FAILURE;
  }

  // (3) Generate heatmaps

  printf("\n(2) GENERATE HEATMAPS\n");
  printf("---------------------\n");

  if (!generate_heatmaps()) {
    perror("Problem in generating heatmaps\n");
    return EXIT_FAILURE;
  }

  // (4) Done

  printf("\nDone.\n");

  return EXIT_SUCCESS;
}

// end
