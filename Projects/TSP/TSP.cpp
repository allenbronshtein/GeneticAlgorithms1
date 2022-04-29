#include "TSP.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include <random>
#include "genlib.h"

static std::vector<std::tuple<double, double>> g_cities;
static size_t g_chromosome_size = 0;
static double g_max_distance = 0;
static double g_min_distance = std::numeric_limits<double>::infinity();
static double g_max_travel_distance = 0;
static double g_min_travel_distance = 0;

void write_to_file(const std::string& file_path, const std::string& content) {
	std::ofstream out(file_path, std::ios::out);
	out << content;
}

template <typename Gen>
void genetic_round(gen_ns::GeneticModel<Gen> genetic_model) {
	std::string file_prefix = GENETIC_OUTPUT_PREFIX;
	std::string file_suffix = GENETIC_OUTPUT_SUFFIX;
	genetic_model.start_sync();
	auto means = genetic_model.get_saved_fitness_mean();
	auto highests = genetic_model.get_saved_highest_fitness();
	auto lowests = genetic_model.get_saved_lowest_fitness();
	size_t size = means.size();
	std::string out = "";
	out =
		std::to_string(genetic_model.get_population_size()) + "\n" +
		std::to_string(genetic_model.get_crossover_prob()) + "\n" +
		std::to_string(genetic_model.get_gen_mutation_prob()) + "\n";
	std::string data = "";
	for (int i = 0; i < size; i++) {
		highests[i] = g_max_travel_distance - highests[i];
		lowests[i] = g_max_travel_distance - lowests[i];
		means[i] = g_max_travel_distance - means[i];
		data += std::to_string(highests[i]) + "," + std::to_string(lowests[i]) + "," + std::to_string(means[i]) + "\n";
	}
	out += data;
	std::string file_name = file_prefix + std::to_string(genetic_model.get_population_size()) + "-" + std::to_string(genetic_model.get_crossover_prob()) + "-" + std::to_string(genetic_model.get_gen_mutation_prob()) + file_suffix;
	write_to_file(file_name, out);
}

void load_cities(const std::string& file_path, std::vector<std::tuple<double, double>>& file_content) {
	std::string fileText;

	std::ifstream file(file_path);
	while (std::getline(file, fileText)) {
		std::vector<double> vect;

		std::stringstream ss(fileText);

		for (int i; ss >> i;) {
			vect.push_back(i);
			if (ss.peek() == ' ')
				ss.ignore();
		}
		file_content.push_back(std::make_tuple(vect[0], vect[1]));
	}

	file.close();
}

double distance(std::tuple<double, double> a, std::tuple<double, double> b) {
	return sqrt(pow(std::get<0>(a) - std::get<0>(b), 2) +
		pow(std::get<1>(a) - std::get<1>(b), 2)) * 1.0;
}

double distance(uint_least8_t city1_i, uint_least8_t city2_i) {
	size_t num_of_cities = g_cities.size();
	double d = std::numeric_limits<double>::infinity();
	if (city1_i < num_of_cities && city2_i < num_of_cities) {
		d = distance(g_cities[city1_i], g_cities[city2_i]);
	}
	return d;
}

double get_max_destination(std::vector<std::tuple<double, double>>  cities)
{
	size_t size = cities.size();
	double max_distance = 0;
	for (size_t i = 0; i < size; i++) {
		for (size_t j = i + 1; j < size; j++) {
			double d = distance(cities[i], cities[j]);
			if (d > max_distance) {
				max_distance = d;
			}
		}
	}
	return max_distance;
}

double get_min_destination(std::vector<std::tuple<double, double>>  cities)
{
	size_t size = cities.size();
	double min_distance = std::numeric_limits<double>::infinity();
	for (size_t i = 0; i < size; i++) {
		for (size_t j = i + 1; j < size; j++) {
			double d = distance(cities[i], cities[j]);
			if (d < min_distance) {
				min_distance = d;
			}
		}
	}
	return min_distance;
}

uint_least8_t get_closest_unvisited_city(std::vector<std::tuple<double, double>> cities_coords, std::vector<bool> is_visited, uint_least8_t current_city) {
	std::vector<double> city_distances(cities_coords.size(), std::numeric_limits<double>::infinity());
	uint_least8_t city_index = -1;
	size_t num_cities = cities_coords.size();
	for (int i = 0; i < num_cities; i++)
	{
		if (!is_visited[i]) {
			city_distances[i] = distance(cities_coords[current_city], cities_coords[i]);
		}
	}
	return min_element(city_distances.begin(), city_distances.end()) - city_distances.begin();
}

//Requires size 8 chromosome, index is row num, value is col num
class TSPInitialize : public gen_ns::InitializeFunction<uint_least8_t > {
public:
	void run(gen_ns::Chromosome<uint_least8_t >* chromosome) const override {
		std::vector<uint_least8_t >values;
		for (size_t i = 0; i < g_chromosome_size; i++) {
			values.push_back(i);
		}
		while (!values.empty()) {
			int i = 0;
			if (values.size() != 1) {
				i = gen_ns::random(0, values.size() - 1);
			}
			chromosome->init_push(values.at(i));
			values.erase(values.begin() + i);
		}
	}
};

class TSPFourPointCrossover : public gen_ns::CrossoverFunction<uint_least8_t> {
public:
	void run(const gen_ns::Chromosome<uint_least8_t>& chromosome1, const gen_ns::Chromosome<uint_least8_t>& chromosome2, std::vector<gen_ns::Chromosome<uint_least8_t>*>& children) const override {
		int first_point = gen_ns::random(0, std::floor(g_chromosome_size / 4) - 1);
		int second_point = gen_ns::random(std::floor(g_chromosome_size / 4) + 1, (std::floor((2*g_chromosome_size )/ 4) - 1));
		int third_point = gen_ns::random((std::floor((2 * g_chromosome_size) / 4) + 1), (std::floor((3 * g_chromosome_size) / 4) + 1));
		int fourth_point = gen_ns::random((std::floor((3 * g_chromosome_size) / 4) + 1),g_chromosome_size-1);
		std::set<uint_least8_t>missing_child1_genes, missing_child2_genes;
		for (size_t i = 0; i < g_chromosome_size; i++) {
			missing_child1_genes.insert(i);
		}
		for (size_t i = 0; i < g_chromosome_size; i++) {
			missing_child2_genes.insert(i);
		}
		std::vector< gen_ns::Chromosome<uint_least8_t>*> _children;
		gen_ns::Chromosome<uint_least8_t>* child1 = new gen_ns::Chromosome<uint_least8_t>(g_chromosome_size);
		gen_ns::Chromosome<uint_least8_t>* child2 = new gen_ns::Chromosome<uint_least8_t>(g_chromosome_size);
		uint_least8_t gen_father, gen_mother;
		double p = gen_ns::prob();
		//start
		for (int i = 0; i <= first_point; i++) {
			if (p < 0.5) {
				gen_father = chromosome1.get_gene_at(i);
				gen_mother = chromosome2.get_gene_at(i);
			}
			else {
				gen_father = chromosome2.get_gene_at(i);
				gen_mother = chromosome1.get_gene_at(i);
			}
			gen_father = chromosome1.get_gene_at(i);
			gen_mother = chromosome2.get_gene_at(i);
			child1->init_push(gen_father);
			missing_child1_genes.erase(gen_father);
			child2->init_push(gen_mother);
			missing_child2_genes.erase(gen_mother);
		}

		//middle 1
		for (int i = first_point + 1; i < second_point; i++) {
			if (p < 0.5) {
				gen_father = chromosome1.get_gene_at(i);
				gen_mother = chromosome2.get_gene_at(i);
			}
			else {
				gen_father = chromosome2.get_gene_at(i);
				gen_mother = chromosome1.get_gene_at(i);
			}
			gen_father = chromosome1.get_gene_at(i);
			gen_mother = chromosome2.get_gene_at(i);
			if (missing_child1_genes.find(gen_mother) != missing_child1_genes.end()) {
				child1->init_push(gen_mother);
				missing_child1_genes.erase(gen_mother);
			}
			if (missing_child2_genes.find(gen_father) != missing_child2_genes.end()) {
				child2->init_push(gen_father);
				missing_child2_genes.erase(gen_father);
			}
		}

		//middle 2
		for (int i = second_point + 1; i < third_point; i++) {
			if (p < 0.5) {
				gen_father = chromosome1.get_gene_at(i);
				gen_mother = chromosome2.get_gene_at(i);
			}
			else {
				gen_father = chromosome2.get_gene_at(i);
				gen_mother = chromosome1.get_gene_at(i);
			}
			gen_father = chromosome1.get_gene_at(i);
			gen_mother = chromosome2.get_gene_at(i);
			if (missing_child1_genes.find(gen_father) != missing_child1_genes.end()) {
				child1->init_push(gen_father);
				missing_child1_genes.erase(gen_father);
			}
			if (missing_child2_genes.find(gen_mother) != missing_child2_genes.end()) {
				child2->init_push(gen_mother);
				missing_child2_genes.erase(gen_mother);
			}
		}

		//middle 3
		for (int i = third_point + 1; i < fourth_point; i++) {
			if (p < 0.5) {
				gen_father = chromosome1.get_gene_at(i);
				gen_mother = chromosome2.get_gene_at(i);
			}
			else {
				gen_father = chromosome2.get_gene_at(i);
				gen_mother = chromosome1.get_gene_at(i);
			}
			gen_father = chromosome1.get_gene_at(i);
			gen_mother = chromosome2.get_gene_at(i);
			if (missing_child1_genes.find(gen_mother) != missing_child1_genes.end()) {
				child1->init_push(gen_mother);
				missing_child1_genes.erase(gen_mother);
			}
			if (missing_child2_genes.find(gen_father) != missing_child2_genes.end()) {
				child2->init_push(gen_father);
				missing_child2_genes.erase(gen_father);
			}
		}

		//end
		for (int i = third_point; i < g_chromosome_size; i++) {
			if (p < 0.5) {
				gen_father = chromosome1.get_gene_at(i);
				gen_mother = chromosome2.get_gene_at(i);
			}
			else {
				gen_father = chromosome2.get_gene_at(i);
				gen_mother = chromosome1.get_gene_at(i);
			}
			gen_father = chromosome1.get_gene_at(i);
			gen_mother = chromosome2.get_gene_at(i);
			if (missing_child1_genes.find(gen_father) != missing_child1_genes.end()) {
				child1->init_push(gen_father);
				missing_child1_genes.erase(gen_father);
			}
			if (missing_child2_genes.find(gen_mother) != missing_child2_genes.end()) {
				child2->init_push(gen_mother);
				missing_child2_genes.erase(gen_mother);
			}
		}

		//child1 leftovers
		while (!missing_child1_genes.empty()) {
			auto it = missing_child1_genes.begin();
			int i = 0;
			if (missing_child1_genes.size() != 1) {
				i = gen_ns::random(0, missing_child1_genes.size() - 1);
			}
			std::advance(it, i);
			child1->init_push(*it);
			missing_child1_genes.erase(*it);
		}

		//child2 leftovers
		while (!missing_child2_genes.empty()) {
			auto it = missing_child2_genes.begin();
			int i = 0;
			if (missing_child2_genes.size() != 1) {
				i = gen_ns::random(0, missing_child2_genes.size() - 1);
			}
			std::advance(it, i);
			child2->init_push(*it);
			missing_child2_genes.erase(*it);
		}

		children.push_back(child1);
		children.push_back(child2);
	}
};

class TSPRouletteWheelSelection : public gen_ns::SelectionFunction<uint_least8_t> {
public:
	void run(const std::vector<gen_ns::Chromosome<uint_least8_t>*>& population, std::queue< int>& parents_index) override{
		size_t population_size = population.size();
		double sum = 0, r;
		for (unsigned int i = 0; i < population_size; i++) {
			sum += population[i]->get_fitness();
		}
		std::uniform_real_distribution<double> unif(0, sum);
		std::default_random_engine re;
		double number_of_selections = std::ceil(population_size * SelectionFunction<uint_least8_t>::pm_select_ratio) - 1;
		for (unsigned k = 0; k < number_of_selections; k++) {
			for (unsigned int j = 0; j < 2; j++) {
				r = unif(re);
				double partial_sum = 0;
				for (unsigned int i = 0; i < population_size; i++) {
					partial_sum += population[i]->get_fitness();
					if (partial_sum >= r) {
						parents_index.push(i);
						break;
					}
				}
			}
		}
	}
};

//Fitness Min Travel distance - Max travel distance
class TSPFitness : public gen_ns::FitnessFunction<uint_least8_t> {
public:
	void run(gen_ns::Chromosome<uint_least8_t>* chromsome) const override {
		double travel_distance = 0;
		for (unsigned int i = 0; i < g_chromosome_size; i++) {
			if (i == g_chromosome_size - 1) {
				travel_distance += distance(chromsome->get_gene_at(i), chromsome->get_gene_at(0));
			}
			else {
				travel_distance += distance(chromsome->get_gene_at(i), chromsome->get_gene_at(i + 1));
			}
		}
		chromsome->set_fitness(g_max_travel_distance - travel_distance);
	}
};

//Mutation swaps values with random index
class TSPMutation : public gen_ns::MutationFunction<uint_least8_t> {
public:
	void run(gen_ns::Chromosome<uint_least8_t>* chromosome) const override {
		double mutation_prob = get_mutation_prob();
		for (unsigned int i = 0; i < g_chromosome_size; i++) {
			if (gen_ns::prob() <= mutation_prob) {
				uint_least8_t temp = chromosome->get_gene_at(i);
				unsigned int j = gen_ns::random(0, g_chromosome_size-1);
				chromosome->set_gene_at(i, chromosome->get_gene_at(j));
				chromosome->set_gene_at(j, temp);
			}
		}
	}
};

void greedy_algorithm() {
	srand(time(NULL));
	int total_cost = 0;
	int number_of_cities = g_cities.size();
	std::vector<bool> is_visited(number_of_cities, false);
	int current_city = rand() % number_of_cities;
	int start_city = current_city;
	is_visited[current_city] = true;
	std::string out = "";
	for (int i = 0; i < number_of_cities - 1; i++)
	{
		int closest_city = get_closest_unvisited_city(g_cities, is_visited, current_city);
		total_cost += distance(g_cities[current_city], g_cities[closest_city]);
		is_visited[closest_city] = true;
		out+= std::to_string(closest_city + 1) + "\n";
		current_city = closest_city;
	}
	out += std::to_string(start_city + 1);
	write_to_file(GREEDY_PATH_OUTPUT_FILE, out);
	total_cost += distance(g_cities[current_city], g_cities[start_city]);
	write_to_file(GREEDY_LENGTH_OUTPUT_FILE, std::to_string(total_cost));
}


int main() {
	load_cities(COORDINATES_FILE,g_cities);
	size_t num_of_cities = g_cities.size();
	g_chromosome_size = num_of_cities;
	g_max_distance = get_max_destination(g_cities);
	g_min_distance = get_min_destination(g_cities);
	g_max_travel_distance = g_max_distance * num_of_cities;
	g_min_travel_distance = g_min_distance * num_of_cities;
	std::cout << "Now executing greedy algo" << std::endl;
	greedy_algorithm();
	std::cout << "Now executing Genetic algo" << std::endl;
	TSPInitialize initialize_function;
	TSPFourPointCrossover crossover_function;
	TSPFitness fitness_function;
	TSPRouletteWheelSelection selection_function;
	TSPMutation mutation_function;
	for (unsigned int population_size = 128; population_size <= 1024; population_size*=2) {
		for (double crossover_prob = 0.5; crossover_prob <= 1; crossover_prob+=0.1) {
			for (double mutation_prob = 0.01; mutation_prob <= 1.2; mutation_prob+=0.5) {
				gen_ns::GeneticModel<uint_least8_t>genetic_model(g_chromosome_size, population_size, &initialize_function, &selection_function, &crossover_function, &fitness_function, &mutation_function, 500000);
				genetic_model.set_generation_limit(1000);
				genetic_model.set_save_generations(20);
				genetic_model.set_crossover_prob(crossover_prob);
				genetic_model.set_gen_mutation_prob(mutation_prob);
				genetic_round(genetic_model);
			}
		}
	}
}