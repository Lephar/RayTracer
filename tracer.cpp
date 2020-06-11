//clang++ tracer.cpp -o tracer -lpthread -w -O2

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <glm/vec3.hpp>
#include <glm/geometric.hpp>

#define EPSILON 0.0009765625

struct ray {
	glm::dvec3 orig;
	glm::dvec3 dir;
} typedef Ray;

struct light {
	glm::dvec3 pos;
	glm::dvec3 col;
	glm::dvec3 aten;
} typedef Light;

struct texture {
	double amb;
	double dif;
	double spec;
	double shin;
	double refl;
} typedef Texture;

struct sphere {
	int pgm;
	int tex;
	glm::dvec3 cent;
	double rad;
} typedef Sphere;

glm::dvec3 traceRay(Ray, int);
glm::dvec3 traceLight(Ray, int, double, int);

pthread_t *threads;
char *image;
int width, height, depth, threadCount, lightCount, pigmentCount, textureCount, sphereCount;
double aspect, fovy, w, h;
glm::dvec3 eye, at, up, fwd, left, *pigments, *pixels;
Light *lights;
Texture *textures;
Sphere *spheres;
unsigned char *data;

void initialize(int argc, char **argv)
{
	if(argc < 2)
	{
		printf("Usage: %s <filename>\n", argv[0]);
		exit(EXIT_SUCCESS);
	}

	if(argc < 3)
		depth = 4;
	else
		depth = atoi(argv[2]);

	if(argc < 4)
		threadCount = get_nprocs();
	else
		threadCount = atoi(argv[3]);
	threads = (pthread_t*) malloc(threadCount * sizeof(pthread_t));

	FILE *input = fopen(argv[1], "r");

	if(input == NULL)
	{
		printf("Cannot open specified input file\n");
		exit(EXIT_FAILURE);
	}

	char placeholder[256];
	image = (char*) malloc(256);

	fscanf(input, "%s", image);
	fscanf(input, "%d %d", &width, &height);
	aspect = (double) width / (double) height;
	pixels = (glm::dvec3*) malloc(height * width * sizeof(glm::dvec3));

	fscanf(input, "%lf %lf %lf", &eye.x, &eye.y, &eye.z);
	fscanf(input, "%lf %lf %lf", &at.x, &at.y, &at.z);
	fscanf(input, "%lf %lf %lf", &up.x, &up.y, &up.z);
	fscanf(input, "%lf", &fovy);

	up = glm::normalize(up);
	fwd = glm::normalize(at - eye);
	left = glm::normalize(glm::cross(up, fwd));
	up = glm::normalize(glm::cross(fwd, left));
	at = eye + fwd;
	h = 2 * tanf(fovy / 2);
	w = aspect * h;

	fscanf(input, "%d", &lightCount);
	lights = (Light*) malloc(lightCount * sizeof(Light));

	for(int i = 0; i < lightCount; i++)
	{
		fscanf(input, "%lf %lf %lf", &lights[i].pos.x, &lights[i].pos.y, &lights[i].pos.z);
		fscanf(input, "%lf %lf %lf", &lights[i].col.r, &lights[i].col.g, &lights[i].col.b);
		fscanf(input, "%lf %lf %lf", &lights[i].aten[0], &lights[i].aten[1], &lights[i].aten[2]);
	}

	fscanf(input, "%d", &pigmentCount);
	pigments = (glm::dvec3*) malloc(pigmentCount * sizeof(glm::dvec3));

	for(int i = 0; i < pigmentCount; i++)
		fscanf(input, "%s %lf %lf %lf", placeholder, &pigments[i].r, &pigments[i].g, &pigments[i].b);

	fscanf(input, "%d", &textureCount);
	textures = (Texture*) malloc(textureCount * sizeof(Texture));

	for(int i = 0; i < textureCount; i++)
		fscanf(input, "%lf %lf %lf %lf %lf", &textures[i].amb, &textures[i].dif, &textures[i].spec, &textures[i].shin, &textures[i].refl);

	fscanf(input, "%d", &sphereCount);
	spheres = (Sphere*) malloc(sphereCount * sizeof(Sphere));

	for(int i = 0; i < sphereCount; i++)
	{
		fscanf(input, "%d %d %s", &spheres[i].pgm, &spheres[i].tex, placeholder);
		fscanf(input, "%lf %lf %lf", &spheres[i].cent.x, &spheres[i].cent.y, &spheres[i].cent.z);
		fscanf(input, "%lf", &spheres[i].rad);
	}

	fclose(input);
}

int compare(double a, double b)
{
	double d = a - b;

	if(d < -EPSILON)
		return -1;
	else if(d > EPSILON)
		return 1;
	else
		return 0;
}

glm::dvec3 intensity(int index, double distance)
{
	return lights[index].col / (lights[index].aten.x + distance * lights[index].aten.y + distance * distance * lights[index].aten.z);
}

glm::dvec3 pigment(int index)
{
	return pigments[spheres[index].pgm];
}

Texture texture(int index)
{
	return textures[spheres[index].tex];
}

double intersect(Ray ray, Sphere sphere)
{
	glm::dvec3 d = ray.orig - sphere.cent;
	double a = glm::dot(ray.dir, ray.dir);
	double b = 2 * glm::dot(ray.dir, d);
	double c = glm::dot(d, d) - sphere.rad * sphere.rad;
	double dsc = b * b - 4 * a * c;

	int state = compare(dsc, 0.0);
	double result = INFINITY;

	if(state == 0)
	{
		double t = -b / (2 * a);

		if(compare(0.0, t) < 0 && compare(t, INFINITY) < 0)
			result = t;
	}

	else if(state == 1)
	{
		double s = sqrt(dsc);
		double t0 = (-b - s) / (2 * a);
		double t1 = (-b + s) / (2 * a);

		if(compare(0.0, t0) < 0 && compare(t0, INFINITY) < 0)
			result = t0;
		else if(compare(0.0, t1) < 0 && compare(t1, result) < 0)
			result = t1;
	}

	return result;
}

glm::dvec3 traceLight(Ray ray, int index, double t, int step)
{
	glm::dvec3 point(ray.orig + ray.dir * (t - EPSILON));
	glm::dvec3 normal(glm::normalize(point - spheres[index].cent));
	glm::dvec3 color(pigment(index) * texture(index).amb * lights[0].col);

	if(glm::distance(ray.orig, spheres[index].cent) < spheres[index].rad)
		normal = -normal;

	for(int i = 1; i < lightCount; i++)
	{
		double tn = INFINITY;
		glm::dvec3 direction(lights[i].pos - point);
		double distance = glm::length(direction);
		Ray light{point, glm::normalize(direction)};

		for(int j = 0; j < sphereCount; j++)
		{
			double z = intersect(light, spheres[j]);

			if(compare(z, tn) < 0 && compare(z, distance) < 0)
				tn = z;
		}

		if(compare(tn, INFINITY) >= 0)
		{
			glm::dvec3 diffuse(texture(index).dif * pigment(index) * glm::max(glm::dot(normal, light.dir), 0.0));
			glm::dvec3 specular(texture(index).spec * powf(glm::max(glm::dot(glm::normalize(ray.orig - point), glm::reflect(-light.dir, normal)), 0.0), texture(index).shin));
			color += intensity(i, distance) * (diffuse + specular);
		}
	}

	if(compare(0.0, texture(index).refl) >= 0)
		return color;

	if(step <= 0)
		return (1.0 - texture(index).refl) * color;

	ray.orig = point;
	ray.dir = glm::reflect(ray.dir, normal);

	return (1.0 - texture(index).refl) * color + texture(index).refl * traceRay(ray, step - 1);
}

glm::dvec3 traceRay(Ray ray, int step)
{
	int index = -1;
	double t = INFINITY;

	for(int i = 0; i < sphereCount; i++)
	{
		double z = intersect(ray, spheres[i]);

		if(compare(z, t) < 0)
		{
			t = z;
			index = i;
		}
	}

	if(compare(t, INFINITY) >= 0)
		return glm::dvec3(0.5, 0.5, 0.5);

	return traceLight(ray, index, t, step);
}

void *calculate(void *threadNumber)
{
	int step = width * height / threadCount;
	int begin = (long) threadNumber * step, end;

	if((long) threadNumber != threadCount - 1)
		end = begin + step;
	else
		end = width * height;

	for(int i = begin; i < end; i++)
	{
		double x = w * (i % width) / width - w / 2;
		double y = h * (i / width) / height - h / 2;

		glm::dvec3 origin(at - x * left - y * up);
		glm::dvec3 direction(glm::normalize(origin - eye));

		pixels[i] = traceRay(Ray{origin, direction}, depth);
	}

	return NULL;
}

void finalize()
{
	data = (unsigned char*) malloc(height * width * 3);

	for(int i = 0; i < height * width; i++)
	{
		data[i * 3 + 0] = (unsigned char) glm::clamp(pixels[i].r * 255.0, 0.0, 255.0);
		data[i * 3 + 1] = (unsigned char) glm::clamp(pixels[i].g * 255.0, 0.0, 255.0);
		data[i * 3 + 2] = (unsigned char) glm::clamp(pixels[i].b * 255.0, 0.0, 255.0);
	}

	free(lights);
	free(pigments);
	free(textures);
	free(spheres);
	free(pixels);

	FILE *output = fopen(image, "wb");

	if(output == NULL)
	{
		printf("Cannot create specified output file\n");
		exit(EXIT_FAILURE);
	}

	fprintf(output, "P6\n%d %d\n255\n", width, height);
	fwrite(data, 1, width * height * 3, output);
	fclose(output);

	free(threads);
	free(image);
	free(data);
}

int main(int argc, char **argv)
{
	initialize(argc, argv);

	for(long i = 0; i < threadCount; i++)
		pthread_create(&threads[i], NULL, calculate, (void*) i);
	for(long i = 0; i < threadCount; i++)
		pthread_join(threads[i], NULL);

	finalize();
}
