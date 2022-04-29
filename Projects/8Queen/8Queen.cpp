#include "8Queen.h"
#include <vector>
#include "genlib.h"
#include <iostream>
#include <ctime>
#include <fstream>
#include <string>

bool is_solved(int chess_board[]) {
	for (int i = 0; i < CHROMOSOME_SIZE - 1; i++)
	{
		for (int j = i + 1; j < CHROMOSOME_SIZE; j++)
		{
			int first_queen_x = i;
			int first_queen_y = chess_board[i];
			int second_queen_x = j;
			int second_queen_j = chess_board[j];

			if (first_queen_y - first_queen_x == second_queen_j - second_queen_x || first_queen_x + first_queen_y == second_queen_j + second_queen_x) {
				return false;
			}
		}
	}

	return true;
}

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
		highests[i] = MAX_COLLISIONS - highests[i];
		lowests[i] = MAX_COLLISIONS - lowests[i];
		means[i] = MAX_COLLISIONS - means[i];
		data += std::to_string(highests[i]) + "," + std::to_string(lowests[i]) + "," + std::to_string(means[i]) + "\n";
	}
	out += data;
	std::string file_name = file_prefix + std::to_string(genetic_model.get_population_size()) + "-" + std::to_string(genetic_model.get_crossover_prob()) + "-" + std::to_string(genetic_model.get_gen_mutation_prob()) + file_suffix;
	write_to_file(file_name, out);
}

void brute_force_eight_queen() {
	std::string out = "";
	auto start = std::chrono::high_resolution_clock::now();
	int chess_board[] = { 0,1,2,3,4,5,6,7 };
	int solutions_found = 0;
	do {
		if (is_solved(chess_board)) {
			for (int i = 0; i < CHROMOSOME_SIZE; i++) {
				std::cout << chess_board[i] << "  ";
			}
			std::cout << std::endl;
			solutions_found++;
			auto stop = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
			out += "Time taken for solution number " + std::to_string(solutions_found) + " : "  + std::to_string(duration.count()) + " milliseconds\n" ;
		}

	} while (std::next_permutation(chess_board, chess_board + CHROMOSOME_SIZE) || solutions_found < NUM_OF_MAX_SOLUTIONS);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	out += "Finished exuecute in: " + std::to_string(duration.count()) + " milliseconds\n";
	write_to_file(BF_OUTPUT_FILE, out);
}

//Requires size 8 chromosome, index is row num, value is col num
class EightQueenInitialize : public gen_ns::InitializeFunction<uint_least8_t > {
public:
	void run(gen_ns::Chromosome<uint_least8_t >* chromosome) const override {
		std::vector<uint_least8_t >values{ 0,1,2,3,4,5,6,7 };
		while (!values.empty()) {
			int i = 0;
			if (values.size() != 1) {
				i = rand() % (values.size() - 1);
			}
			chromosome->init_push(values.at(i));
			values.erase(values.begin() + i);
		}
	}
};

class EightQueenCrossover : public gen_ns::CrossoverFunction<uint_least8_t> {
public:
	void run(const gen_ns::Chromosome<uint_least8_t>& chromosome1, const gen_ns::Chromosome<uint_least8_t>& chromosome2, std::vector<gen_ns::Chromosome<uint_least8_t>*>& children) const override {
		std::vector<gen_ns::Chromosome<uint_least8_t>> parents;
		parents.push_back(chromosome1);
		parents.push_back(chromosome2);
		std::set<uint_least8_t >missing_child1_genes{ 0,1,2,3,4,5,6,7 };
		std::set<uint_least8_t >missing_child2_genes{ 0,1,2,3,4,5,6,7 };
		std::vector< std::set<uint_least8_t >*> _missing_children_genes;
		_missing_children_genes.push_back(&missing_child1_genes);
		_missing_children_genes.push_back(&missing_child2_genes);
		std::vector< gen_ns::Chromosome<uint_least8_t>*> _children;
		gen_ns::Chromosome<uint_least8_t>* child1 = new gen_ns::Chromosome<uint_least8_t>(CHROMOSOME_SIZE);
		gen_ns::Chromosome<uint_least8_t>* child2 = new gen_ns::Chromosome<uint_least8_t>(CHROMOSOME_SIZE);
		_children.push_back(child1);
		_children.push_back(child2);
		for (unsigned int child_i = 0; child_i < 2; child_i++) { // for each child
			gen_ns::Chromosome< uint_least8_t>* child = _children.at(child_i);
			std::set<uint_least8_t >* missing_genes = _missing_children_genes.at(child_i);

			for (unsigned int gen_i = 0; gen_i < CHROMOSOME_SIZE; gen_i++) { // for each gen
				int parent_i = 0;
				if (gen_ns::prob() > 0.5) {
					parent_i = 1;
				}
				const gen_ns::Chromosome<uint_least8_t> parent = parents.at(parent_i);
				uint_least8_t gen = parent.get_gene_at(gen_i);
				if (missing_genes->find(gen) != missing_genes->end()) { // child doesn't have chosen parent chromosome
					child->init_push(gen);
					missing_genes->erase(gen);
					continue;
				}
				const gen_ns::Chromosome<uint_least8_t>&alt_parent = parents.at((parent_i + 1) % 2);
				gen = alt_parent.get_gene_at(gen_i);
				if (missing_genes->find(gen) != missing_genes->end()) {// child has first chosen parent chromosome, but not second
					child->init_push(gen);
					missing_genes->erase(gen);
					continue;
				}
				else {// child has both parents gen, so we choose randomly
					int j = 0;
					if (missing_genes->size() != 1) {
						j = rand() % (missing_genes->size() - 1);
					}
					std::set<uint_least8_t >::iterator it = missing_genes->begin();
					std::advance(it, j);
					child->init_push(*it);
					missing_genes->erase(*it);
				}
			}// for gen
			children.push_back(child);
		}// for child
	}
};

//Fitness 0-56
class EightQueenFitness : public gen_ns::FitnessFunction<uint_least8_t> {
public:
	void run(gen_ns::Chromosome<uint_least8_t>* chromsome) const override {
		unsigned int num_of_collisions = 0;
		for (unsigned int i = 0; i < CHROMOSOME_SIZE; i++) {
			for (unsigned j = i + 1; j < CHROMOSOME_SIZE; j++) {
				if (i + chromsome->get_gene_at(i) == j + chromsome->get_gene_at(j) || i - chromsome->get_gene_at(i) == j - chromsome->get_gene_at(j)) {
					num_of_collisions += 2;
				}
			}
		}
		chromsome->set_fitness((MAX_COLLISIONS - num_of_collisions));
	}
};

void print_solutions(const std::vector<gen_ns::Chromosome<uint_least8_t>*>& solutions) {
	size_t num_of_solutions = solutions.size();
	for (unsigned int i = 0; i < num_of_solutions; i++) {
		printf("Solution %d : ", i + 1);
		gen_ns::Chromosome<uint_least8_t>* solution = solutions.at(i);
		const std::vector<uint_least8_t>& genes = solution->get_genes();
		for (unsigned int j = 0; j < CHROMOSOME_SIZE; j++) {
			uint_least8_t col = genes.at(j);
			ROW;
			for (unsigned k = 0; k <= CHROMOSOME_SIZE; k++) {
				COL;
				if (col == k) {
					printf("   Q");
				}
			}
		}
		ROW;
	}
}

//Mutation swaps values with random index
class EightQueenMutation : public gen_ns::MutationFunction<uint_least8_t> {
public:
	void run(gen_ns::Chromosome<uint_least8_t>* chromosome) const override {
		double mutation_prob = get_mutation_prob();
		for (unsigned int i = 0; i < CHROMOSOME_SIZE; i++) {
			if (gen_ns::prob() <= mutation_prob) {
				uint_least8_t temp = chromosome->get_gene_at(i);
				unsigned int j = rand() % CHROMOSOME_SIZE;
				chromosome->set_gene_at(i, chromosome->get_gene_at(j));
				chromosome->set_gene_at(j, temp);
			}
		}
	}
};

int main()
{
	system("mkdir .\\results\\brute_force\\");
	system("mkdir .\\results\\genetic\\");
	std::cout << "Now executing brute force" << std::endl;
	brute_force_eight_queen();
	std::cout << "\nNow executing Genetic algo" << std::endl;
	EightQueenInitialize initialize_function;
	EightQueenCrossover crossover_function;
	EightQueenFitness fitness_function;
	gen_ns::selection_func_ns::ExploreExploitSelection< uint_least8_t> selection_function;
	EightQueenMutation mutation_function;

	//modify values
	//for (unsigned int population_size = 128; population_size <= 1024; population_size *= 2) {
	//	for (double crossover_prob = 0.5; crossover_prob <= 1; crossover_prob += 0.1) {
	//		for (double mutation_prob = 0.01; mutation_prob <= 1.2; mutation_prob += 0.5) {

	//------------------------------------------------User values------------------------------------------//
	unsigned int population_size = 1024;
	double crossover_prob = 0.8;
	double mutation_prob = 1;
	unsigned int record_data_per_generations = 20; //Per how much generations to record data
	unsigned int generation_limit = 1000; //After how much generations program terminates
	//----------------------------------------------------------------------------------------------------------//
				gen_ns::GeneticModel<uint_least8_t>genetic_model(CHROMOSOME_SIZE, population_size, &initialize_function, &selection_function, &crossover_function, &fitness_function, &mutation_function, MAX_COLLISIONS);
				genetic_model.set_min_number_of_solutions(92);
				genetic_model.set_generation_limit(generation_limit);
				genetic_model.set_save_generations(record_data_per_generations);
				genetic_model.set_crossover_prob(crossover_prob);
				genetic_model.set_gen_mutation_prob(mutation_prob);
				genetic_round(genetic_model);
	//		}
	//	}
	//}
}