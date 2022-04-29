#ifndef GEN_LIB_H
#define GEN_LIB_H
#include <vector>
#include <iostream>
#include <random>
#include <cmath>
#include <queue>
#include <set>
#include <chrono>

#define DEFAULT_CROSSOVER_PROB 0.7 //Probability of a couple to make children
#define DEFAULT_GEN_MUTATION_PROB 0.01 //Probability of a gene to mutate
#define DEFAULT_CHROMOSOME_MUTATION_PROB 1 //Probability of a chromosome to be up for mutation
#define DEFAULT_EXPLOIT_PROB 1 //Probability strong will mate only with strong
#define DEFAULT_SELECT_RATIO 0.4 // 0 == No couples are created 1 == (population size - 1) couples are created
#define DEFAULT_ERADICATE_SLOPE 100 //Min diff required to prevent eradicate
#define DEFAULT_ERADICATE_GENERATIONS 100 //Number of generations to cooldown between eradicates
#define DEFAULT_ERADICATE_PERCENTAGE 0.1 //How much of population % do we want to eradicate
#define DEFAULT_ERADICATE_PROB 0 //When eradicate generation comes, probability eradicate will actually happen
#define DEFAULT_AFTER_GENERATIONS_PRINT 1 //After this much generations has passed from previous print, we will print again
#define DEFAULT_SAVE_MEAN_AFTER_GENERATIONS 100 //After this much generations has passed from previous print, we will print again

namespace gen_ns
{
	enum {
		OFF,
		ON,
	}typedef run_t;
	enum {
		WEAK,
		STRONG,
	}typedef rank_t;
	/* 
	chooses a random value between min_val and max_val, if value is int max value should be + 1 
	*/
	int random(int min_val, int max_val) {
		std::random_device rd;
		std::mt19937 mt(rd());
		std::uniform_int_distribution<int> dist(min_val, max_val);
		int random_number = dist(mt);
		return random_number;
	}

	//chooses random number between 0 to 1
	static double prob() {
		std::uniform_real_distribution<double> unif(0, 1);
		std::default_random_engine re;
		double a_random_double = unif(re);
		return a_random_double;
	}

	template<typename Gen>
	class Chromosome
	{
		typedef std::vector<Gen> Genes;
	private:
		bool m_is_max_set_ = false, m_is_min_set_ = false, m_is_init_ = false, m_is_fitness_set_ = false;
		size_t m_size_;
		Gen m_max_value_, m_min_value_;
		Genes m_genes_;
		double m_fitness_ = 0;
	public:
		// in constructor we require size of chromosome
		Chromosome(size_t size)
		{
			m_size_ = size;
			m_max_value_ = std::numeric_limits<Gen>::max();
			m_min_value_ = std::numeric_limits<Gen>::min();
		}

		//Clears chromosome
		void clear_genes() {
			m_genes_.clear();
			m_is_init_ = false;
			m_fitness_ = 0;
		}

		//Adds new genes to chromosome. Can only be done if m_is_init_ is false, when chromosome is full m_is_init_ turns to true
		void init_push(Gen gen) {
			size_t num_genes = m_genes_.size();
			if (!m_is_init_&&num_genes < m_size_) {
				m_genes_.push_back(gen);
			}
			num_genes = m_genes_.size();
			if (num_genes == m_size_) {
				m_is_init_ = true;
			}
		}

		//returns is chromosome init
		bool is_init() const {
			return m_is_init_;
		}

		//sets is init. can only set once from false to true, once true it cannot change
		void set_is_init() {
			m_is_init_ = true;
		}

		// returns size 
		size_t get_size() const
		{
			return m_size_;
		}

		// returns max value
		const Gen& get_max_value() const
		{
			return m_max_value_;
		}

		//returns min value
		const Gen& get_min_value() const
		{
			return m_min_value_;
		}

		// returns vector of genes
		const Genes& get_genes() const
		{
			return m_genes_;
		}

		//returns gen at specified index, if index is bad returns null
		const Gen& get_gene_at(unsigned int index) const
		{
			Gen* value = nullptr;
			if (index < m_size_) {
				return m_genes_.at(index);
			}
			return *value;
		}

		//set the max value available. can be set only once, changes genes accordingly
		void set_max_value(const Gen& max_value)
		{
			if (!m_is_max_set_) {
				m_max_value_ = max_value;
				if (m_is_init_) {
					for (unsigned int i = 0; i < m_size_; i++) {
						if (m_genes_.at(i) > m_max_value_) {
							m_genes_.at(i) = m_max_value_;
						}
					}
				}
				m_is_max_set_ = true;
			}
		}

		//set the min value available. can be set only once, changes genes accordingly
		void set_min_value(const Gen& min_value)
		{
			if (!m_is_min_set_) {
				m_min_value_ = min_value;
				if (m_is_init_) {
					for (unsigned int i = 0; i < m_size_; i++) {
						if (m_genes_.at(i) < m_min_value_) {
							m_genes_.at(i) = m_min_value_;
						}
					}
				}
				m_is_min_set_ = true;
			}
		}

		//sets all genes  to given value
		void set_genes(const Gen& value)
		{
			for (unsigned int i = 0; i < m_size_; i++) {
				if (value <= m_max_value_ && m_min_value_ <= value) {
					m_genes_.at(i) = value;
				}
			}
		}

		//sets a specified gen to a given value
		void set_gene_at(unsigned int index, const Gen& value)
		{
			if (index < m_size_ && value <= m_max_value_ && m_min_value_ <= value) {
				m_genes_.at(index) = value;
			}
		}

		//check if a gen is valid (between max and min), returns false if bad index given
		bool is_gene_valid(unsigned int index) const
		{
			if (index < m_size_) {
				return m_genes_.at(index) <= m_max_value_ && m_genes_.at(index) >= m_min_value_;
			}
			return false;
		}

		//check if current chromosome is valid 
		bool is_chromosome_valid() const
		{
			for (Gen gene : m_genes_) {
				if (!is_gene_valid(gene)) {
					return false;
				}
			}
			return true;
		}

		//returns fitness of chromsome (it doesn't calculated, it returns saved value)
		double get_fitness() const {
			return m_fitness_;
		}

		//sets fitness function of chromsome
		void set_fitness(double fitness) {
			if (fitness < 0) {
				m_fitness_ = 0;
			}
			else {
				m_fitness_ = fitness;
			}
			m_is_fitness_set_ = true;
		}

		//checks if one chromosome is same as other
		bool operator==(const Chromosome & other) const
		{
			if (this->m_max_value_ != other.m_max_value_ || this->m_min_value_ != other.m_min_value_) {
				return false;
			}
			if (m_size_ == other.get_size()) {
				for (unsigned int i = 0; i < m_size_; i++) {
					if (m_genes_.at(i) != other.get_gene_at(i)) {
						return false;
					}
				}
				return true;
			}
			return false;
		}

		//check if one chromsome is different from another
		bool operator!=(const Chromosome & other) const
		{
			return !(this == other);
		}

		//assignment operator
		Chromosome<Gen> & operator=(const Chromosome & other)
		{
			m_is_max_set_ = other.m_is_max_set_;
			m_is_min_set_ = other.m_is_min_set_;
			m_is_init_ = other.m_is_init_;
			m_size_ = other.m_size_;
			m_max_value_ = other.m_max_value_;
			m_min_value_ = other.m_min_value_;
			m_genes_ = other.m_genes_;
			m_fitness_ = other.m_fitness_;
			return *this;
		}

		bool operator<(const Chromosome & other) const
		{
			return this->get_fitness() < other.get_fitness();
		}

		~Chromosome() {	}
	};

	//Chromosome initializer
	template<typename Gen>
	class InitializeFunction {
	public:
		virtual void run(Chromosome<Gen>* chromosome) const = 0;
	};

	//Selection base class, it assumes indexes are sorted from strongest to weakest
	template <typename Gen>
	class SelectionFunction {
	protected:
		double pm_exploit_prob = DEFAULT_EXPLOIT_PROB;
		double pm_select_ratio = DEFAULT_SELECT_RATIO;
	public:
		virtual void run(const std::vector<Chromosome<Gen>*>& population_size, std::queue<int>& population_index_selected) = 0;
		virtual double get_exploit_prob() {
			return pm_exploit_prob;
		}
		virtual double get_select_ratio() {
			return pm_select_ratio;
		}
		virtual void set_exploit_prob(double exploit_prob) {
			if (exploit_prob > 1) {
				pm_exploit_prob = 1;
			}
			else if (exploit_prob < 0) {
				pm_exploit_prob = 0;
			}
			else {
				pm_exploit_prob = exploit_prob;
			}
		}
		virtual void set_select_ratio(double select_ratio) {
			if (select_ratio > 1) {
				pm_select_ratio = 1;
			}
			else if (select_ratio < 0) {
				pm_select_ratio = 0;
			}
			else {
				pm_select_ratio = select_ratio;
			}
		}
	};

	//Crossover base class
	template<typename Gen>
	class CrossoverFunction {
	protected:
		double pm_crossover_prob = DEFAULT_CROSSOVER_PROB;
	public:
		virtual double get_crossover_prob() {
			return pm_crossover_prob;
		}
		virtual void set_crossover_prob(double crossover_prob) {
			if (crossover_prob > 1) {
				pm_crossover_prob = 1;
			}
			else if (crossover_prob < 0) {
				pm_crossover_prob = 0;
			}
			else {
				pm_crossover_prob = crossover_prob;
			}
		}
		virtual void run(const Chromosome<Gen>& chromsome1, const Chromosome<Gen>& chromsome2, std::vector<Chromosome<Gen>*>& children) const = 0;
	};

	//Fitness base class
	template<typename Gen>
	class FitnessFunction {
	public:
		virtual void run(Chromosome<Gen>* chromsome) const = 0;
	};

	//Mutation base class
	template<typename Gen>
	class MutationFunction {
	protected:
		double pm_mutation_prob = DEFAULT_GEN_MUTATION_PROB;
	public:
		virtual void run(Chromosome<Gen>* chromosome) const = 0;
		virtual double get_mutation_prob() const {
			return pm_mutation_prob;
		}
		virtual void set_mutation_prob(double mutation_prob) {
			if (mutation_prob > 1) {
				pm_mutation_prob = 1;
			}
			else if (mutation_prob < 0) {
				pm_mutation_prob = 0;
			}
			else {
				pm_mutation_prob = mutation_prob;
			}
		}
	};

	template<typename Gen>
	class GeneticModel {
		typedef std::vector<Chromosome<Gen>*> Population;
	private:
		bool m_is_population_init_ = false;
		bool m_is_population_sorted_ = true;

		size_t m_min_number_of_solutions_ = 1;
		unsigned int m_number_of_found_solutions_ = 0;
		Population m_solutions_;
		double m_min_solution_fitness_;
		bool m_found_min_required_solutions_ = false;

		run_t m_is_execution_on_ = run_t::ON;
		bool m_run_ = false;
		const size_t m_chromosome_size_;
		double m_crossover_prob_ = DEFAULT_CROSSOVER_PROB;
		double m_chromosome_mutation_prob_ = DEFAULT_CHROMOSOME_MUTATION_PROB;
		double m_gen_mutation_prob_ = DEFAULT_GEN_MUTATION_PROB;
		Gen m_max_value_;
		Gen m_min_value_;
		bool m_is_max_value_set_ = false;
		bool m_is_min_value_set_ = false;
		double m_exploit_prob_ = DEFAULT_EXPLOIT_PROB;
		double m_select_ratio_ = DEFAULT_SELECT_RATIO;
		size_t m_population_size_;
		Population m_population_;
		InitializeFunction<Gen>* m_initialize_function_;
		SelectionFunction<Gen>* m_selection_function_;
		CrossoverFunction<Gen>* m_crossover_function_;
		FitnessFunction<Gen>* m_fitness_function_;
		MutationFunction<Gen>* m_mutation_function_;
		double m_mean_ = NULL;
		double m_median_ = NULL;
		double m_mean_diff_ = 0;
		unsigned int m_generations_ = 0;
		double m_eradicate_generations_ = DEFAULT_ERADICATE_GENERATIONS;
		double m_eradicate_slope_ = DEFAULT_ERADICATE_SLOPE;
		double m_eradicate_percentage_ = DEFAULT_ERADICATE_PERCENTAGE;
		size_t m_num_of_eradications_ = 0;
		size_t m_previous_recorded_generation_ = 0;
		double m_previous_recorded_mean_ = NULL;
		double m_eradicate_prob_ = DEFAULT_ERADICATE_PROB;
		
		unsigned int m_previous_print_generation_ = 0;
		unsigned int m_print_after_passed_generations_ = DEFAULT_AFTER_GENERATIONS_PRINT;

		unsigned int m_previous_saved_mean_generation_ = 0;
		unsigned int m_save_mean_after_passed_generations_ = DEFAULT_SAVE_MEAN_AFTER_GENERATIONS;

		std::vector<unsigned int>m_saved_generation_;
		std::vector<double>m_saved_fitness_mean_;
		std::vector<double>m_saved_highest_fitness_;
		std::vector<double>m_saved_lowest_fitness_;
		run_t m_execute_parents_mode_ = run_t::OFF;
		run_t m_add_existing_chromosomes_mode_ = run_t::ON;
		unsigned int m_generation_limit_ = std::numeric_limits<unsigned int>::infinity();

	public:
		// in constructor we initialize chromosome size, population size, crossover function, fitness function and mutation function
		GeneticModel(size_t chromosome_size, size_t population_size, InitializeFunction<Gen>* initialize_function, SelectionFunction<Gen>*  selection_function, CrossoverFunction<Gen>*  crossover_function, FitnessFunction<Gen>*  fitness_function, MutationFunction<Gen>* mutation_function, double min_solution_fitness)
			: m_chromosome_size_(chromosome_size)
		{
			if (population_size == 0) {
				m_population_size_ = 1;
			}
			m_population_size_ = population_size;
			m_initialize_function_ = initialize_function;
			m_selection_function_ = selection_function;
			m_crossover_function_ = crossover_function;
			m_fitness_function_ = fitness_function;
			m_mutation_function_ = mutation_function;
			m_max_value_ = std::numeric_limits<Gen>::max();
			m_min_value_ = std::numeric_limits<Gen>::min();
			m_gen_mutation_prob_ = mutation_function->get_mutation_prob();
			m_min_solution_fitness_ = min_solution_fitness;
		}

		//get generations limitiation
		unsigned int get_generation_limit() {
			return m_generation_limit_;
		}

		//sets limit on number of generations to run
		void set_generation_limit(unsigned int limit) {
			m_generation_limit_ = limit;
		}

		//returns highest recorded fitness
		const std::vector<double>& get_saved_highest_fitness() {
			return m_saved_highest_fitness_;
		}

		//returns highest recorded fitness
		const std::vector<double>& get_saved_lowest_fitness() {
			return m_saved_lowest_fitness_;
		}

		//return mode for adding existing chromosome
		run_t get_add_existing_chromosomes_mode() {
			return m_add_existing_chromosomes_mode_;
		}

		//sets mode for adding existing chromosome
		void set_add_existing_chromosomes_mode(run_t mode) {
			m_add_existing_chromosomes_mode_ = mode;
		}
		
		//check if given chromosome is in population
		bool is_chromosome_in_population(const Chromosome<Gen>* chromosome) {
			size_t size = m_population_.size();
			for (unsigned int i = 0; i < size; i++) {
				if (*chromosome == *(m_population_[i])) {
					return true;
				}
			}
			return false;
		}

		//returns execute parents mode
		run_t get_execute_parents_mode() {
			return m_execute_parents_mode_;
		}

		//sets execute parents mode
		void set_execute_parents_mode(run_t mode) {
			m_execute_parents_mode_ = mode;
		}

		//returns eradicate percentage
		double get_eradicate_percentage() {
			return m_eradicate_percentage_;
		}

		//sets eradicate percentage
		void set_eradicate_percentage(double eradicate_percentage) {
			if (eradicate_percentage < 0) {
				m_eradicate_percentage_ = 0;
			}
			else if (eradicate_percentage > 1) {
				m_eradicate_percentage_ = 1;
			}
			else {
				m_eradicate_percentage_ = eradicate_percentage;
			}
		}

		//returns eradicate probability
		double get_eradicate_prob() {
			return m_eradicate_prob_;
		}

		//sets eradicate probability
		void set_eradicate_prob(double eradicate_prob) {
			if (eradicate_prob > 1) {
				m_eradicate_prob_ = 1;
			}
			else if (eradicate_prob < 0) {
				m_eradicate_prob_ = 0;
			}
			else {
				m_eradicate_prob_ = eradicate_prob;
			}
		}

		//returns saved generations
		const std::vector<unsigned int>& get_saved_generations() {
			return m_saved_generation_;
		}

		//returns saved fitness means
		const std::vector<double>& get_saved_fitness_mean() {
			return m_saved_fitness_mean_;
		}

		//return minimum number of solutions
		size_t get_min_number_of_solutions() const {
			return m_min_number_of_solutions_;
		}

		//sets minimum number of solutions. If it's larger than population size, it will be set to population size
		void set_min_number_of_solutions(size_t min_number_of_solutions) {
			if (min_number_of_solutions > m_population_size_) {
				m_min_number_of_solutions_ = m_population_size_;
			}
			else {
				m_min_number_of_solutions_ = min_number_of_solutions;
			}
		}

		//sets after how much generations we print data
		void set_print_generations(unsigned int print_generations) {
			m_print_after_passed_generations_ = print_generations;
		}

		//sets after how much generations we save data
		void set_save_generations(unsigned int save_mean_generations) {
			m_save_mean_after_passed_generations_ = save_mean_generations;
		}

		//returns number of solutions found
		size_t get_number_of_found_solutions() const {
			return m_number_of_found_solutions_;
		}

		//returns min fitness
		double get_min_solution_fitness() const {
			return m_min_solution_fitness_;
		}

		//set min solutions fitness
		void set_min_solution_fitness(double min_fitness) {
			m_min_solution_fitness_ = min_fitness;
		}

		//returns solutions found
		const Population& get_solutions() const {
			return m_solutions_;
		}

		//returns if minmum number of solutions is found
		bool is_min_number_of_solutions_found() {
			return m_found_min_required_solutions_;
		}

		//return exploit probability
		double get_exploit_prob() const {
			return m_exploit_prob_;
		}

		//sets exploit probabilty
		void set_exploit_prob(double exploit_prob) {
			if (exploit_prob > 1) {
				m_exploit_prob_ = 1;
			}
			else if (exploit_prob < 0) {
				m_exploit_prob_ = 0;
			}
			else {
				m_exploit_prob_ = exploit_prob;
			}
			m_selection_function_->set_exploit_prob(exploit_prob);
		}

		//return select ratio
		double get_select_ratio() const {
			return m_select_ratio_;
		}

		//sets select ratio
		void set_select_ratio(double select_ratio)  {
			if (select_ratio > 1) {
				m_select_ratio_ = 1;
			}
			else if (select_ratio < 0) {
				m_select_ratio_ = 0;
			}
			else {
				m_select_ratio_ = select_ratio;
			}
			m_selection_function_->set_select_ratio(select_ratio);
		}

		//returns max size of population
		size_t get_population_size() const {
			return m_population_size_;
		}

		//returns current size of population
		size_t get_current_population_size() const {
			return m_population_.size();
		}

		//sets size of population
		void set_population_size(size_t population_size) {
			m_population_size_ = population_size;
		}

		//returns chromosomes in model
		const Population& get_population() const {
			return m_population_;
		}

		//runs model
		void start_sync() {
			if (!m_is_population_init_) {
				init_population_();
			}
			m_run_ = true;
			thread_function_();
		}

		//stops model
		void stop() {
			size_t over_population_count = 0;
			m_run_ = false;
		}

		//returns probability of mutation
		double get_gen_mutation_prob() const {
			return m_gen_mutation_prob_;
		}

		//sets probability of mutation
		void set_gen_mutation_prob(double mutation_prob) {
			if (mutation_prob > 1) {
				m_gen_mutation_prob_ = 1;
			}
			else if (mutation_prob < 0) {
				m_gen_mutation_prob_ = 0;
			}
			else {
				m_gen_mutation_prob_ = mutation_prob;
			}
			m_mutation_function_->set_mutation_prob(mutation_prob);
		}

		//returns probability of crossover
		double get_crossover_prob() const {
			return m_crossover_prob_;
		}

		//sets probability of crossover
		void set_crossover_prob(double crossover_prob) {
			if (crossover_prob > 1) {
				m_crossover_prob_ = 1;
			}
			else if (crossover_prob < 0) {
				m_crossover_prob_ = 0;
			}
			else {
				m_crossover_prob_ = crossover_prob;
			}
			m_crossover_function_->set_crossover_prob(crossover_prob);
		}

		//returns max value of chromosomes genes
		const Gen& get_max_value() const {
			return m_max_value_;
		}

		//returns min value of chromosomes genes
		const Gen& get_min_value() const {
			return m_min_value_;
		}

		//sets maximum value of chromosomes genes, sets it to all chromosomes in population as well. One time function
		void set_max_value(const Gen& max_value) {
			if (!m_is_max_value_set_) {
				m_max_value_ = max_value;
				if (m_is_population_init_) {
					for (Chromosome<Gen>* chromsome : m_population_) {
						chromsome->set_max_value(max_value);
					}
				}
				m_is_max_value_set_ = true;
			}
		}

		//sets minimum value of chromosomes genes, sets it too all chromosomes in population as well. One time function
		void set_min_value(const Gen& min_value) {
			if (!m_is_min_value_set_) {
				m_min_value_ = min_value;
				if (m_is_population_init_) {
					for (Chromosome<Gen>* chromsome : m_population_) {
						chromsome->set_min_value(min_value);
					}
				}
				m_is_min_value_set_ = true;
			}
		}

		//returns true if over populated. Saves how much in over_population_count. If number is negative it means it's underpopulated
		bool is_over_populated(size_t& over_population_count) const {
			over_population_count = 0;
			if (m_population_.size() > m_population_size_) {
				over_population_count = m_population_.size() - m_population_size_;
				return true;
			}
			return false;
		}

		//returns true if under populated. Saves how much in under_population_count. If number is negative it means it's overpopulated
		bool is_under_populated(size_t& under_population_count) const {
			under_population_count = m_population_size_ - m_population_.size();
			return under_population_count > 0;
		}

		void MergeSortedIntervals(unsigned int start, unsigned int middle, unsigned int end) {
			std::vector<Chromosome<Gen>*> temp;
			int i, j;
			i = start;
			j = middle + 1;
			while (i <= middle && j <= end) {
				if (m_population_.at(i)->get_fitness() >= m_population_.at(j)->get_fitness()) {
					temp.push_back(m_population_.at(i));
					++i;
				}
				else {
					temp.push_back(m_population_.at(j));
					++j;
				}
			}
			while (i <= middle) {
				temp.push_back(m_population_.at(i));
				++i;
			}
			while (j <= end) {
				temp.push_back(m_population_.at(j));
				++j;
			}
			for (int i = start; i <= end; ++i) {
				m_population_.at(i) = temp.at(i - start);
			}
			printf("done");
		}
		//sorts population by fitness [largest, ... ,smallest], if already sorted it will not sort
		void sort_population() {
			if (!m_is_population_sorted_) {
				std::sort(m_population_.begin(), m_population_.end(), [](const Chromosome<Gen>* chromosome1, const Chromosome<Gen>* chromosome2)
				{
					return chromosome1->get_fitness() > chromosome2->get_fitness();
				});
			}
			m_is_population_sorted_ = true;
		}

		//add to population
		void insert_to_population(Chromosome<Gen>* chromosome) {
			m_population_.push_back(chromosome);
			m_is_population_sorted_ = false;
		}

		//add to population sorted
		void insert_to_population_sorted(Chromosome<Gen>* chromosome) {
			if (!m_is_population_sorted_) {
				insert_to_population(chromosome);
				return;
			}
			size_t population_size = m_population_.size();
			for (size_t i = population_size - 1; i >= 0; i--) {
				if (chromosome->get_fitness() <= m_population_.at(i)->get_fitness()) {
					m_population_.insert(m_population_.begin() + i + 1, chromosome);
					return;
				}
			}
			m_population_.insert(m_population_.begin(), chromosome);
		}

		//remove k-th largest from population [1 - current population size]
		void remove_weakest() {
			sort_population();
			delete *(m_population_.end() - 1);
			m_population_.erase(m_population_.end() - 1);
		}

		//removes strongest from population
		void remove_strongest() {
			sort_population();
			delete *(m_population_.end() - 1);
			m_population_.erase(m_population_.begin());
		}

		//remove by index from population
		void remove_by_index(unsigned int i) {
			if (i < m_population_.size()) {
				delete *(m_population_.begin() + i);
				m_population_.erase(m_population_.begin() + i);
			}
		}

		// get k-th largest from population [1- current population size]
		const Chromosome<Gen>* get_kth_chromosome(unsigned int k) {
			sort_population();
			if (k <= m_population_.size()) {
				return m_population_.at(k - 1);
			}
			return NULL;
		}

		//removes k chromosomes
		void execute(size_t k, rank_t rank) {
			if (k > m_population_.size()) {
				k = m_population_.size();
			}
			for (unsigned int i = 0; i < k; i++) {
				if (rank == rank_t::WEAK) {
					remove_weakest();
				}
				else if (rank == rank_t::STRONG) {
					remove_weakest();
				}
			}
			if (m_population_.size() == 0) {
				m_is_population_init_ = false;
				m_is_population_sorted_ = true;
			}
			else if (m_population_.size() == 1) {
				m_is_population_sorted_ = true;
			}
		}

		//removes all chromosomes
		void eradicate(double eradicate_percentage) {
			if (eradicate_percentage > 1) {
				eradicate_percentage = 1;
			}
			else if (eradicate_percentage < 0) {
				eradicate_percentage = 0;
			}
			size_t population_size = m_population_.size();
			double num_to_eradicate = std::floor(population_size*eradicate_percentage);
			if (num_to_eradicate == population_size) {
				m_is_population_init_ = false;
			}
			execute(num_to_eradicate,rank_t::STRONG);
		}

		//returns is execution on
		bool is_execution_on() const {
			return m_is_execution_on_;
		}

		//sets execution on/off
		void set_execution_mode(run_t mode) {
			m_is_execution_on_ = mode;
		}

		//checks if solution exists
		bool solution_exist(Chromosome<Gen>* new_solution) {
			size_t size = m_solutions_.size();
			for (unsigned int i = 0; i < size; i++) {
				if (*(m_solutions_.at(i)) == *new_solution) {
					return true;
				}
			}
			return false;
		}

		//get chromosome mutation probability
		double get_chromosome_mutation_prob() {
			return m_chromosome_mutation_prob_;
		}

		//sets chromosome mutation probability
		void set_chromosome_mutation_prob(double chromosome_mutation_prob) {
			if (chromosome_mutation_prob > 1) {
				m_chromosome_mutation_prob_ = 1;
			}
			else if (chromosome_mutation_prob < 0) {
				m_chromosome_mutation_prob_ = 0;
			}
			else {
				m_chromosome_mutation_prob_ = chromosome_mutation_prob;
			}
		}

		//returns generations passed
		unsigned int get_generations_passed() {
			return m_generations_;
		}

		//returns after how much generations we may eradicate
		double get_eradicate_generations() {
			return m_eradicate_generations_;
		}

		//sets after how much generations we may eradicate
		void set_eradicate_generations(double eradicate_generations) {
			if (eradicate_generations <= 0) {
				m_eradicate_generations_ = 0;
			}
			else {
				m_eradicate_generations_ = std::ceil(eradicate_generations);
			}
		}

		//returns eradicate slope
		double get_eradicate_slope() {
			return m_eradicate_slope_;
		}

		//sets eradicate slope
		void set_eradicate_slope(double eradicate_slope) {
			m_eradicate_slope_ = eradicate_slope;
		}

		//returns number of eradications
		size_t get_num_of_eradications() {
			return m_num_of_eradications_;
		}

		~GeneticModel() {	}

	private:

		// initializes population
		void init_population_() {
			if (m_is_population_init_) {
				return;
			}
			size_t current_population_size = m_population_.size();
			for (size_t i = current_population_size; i < m_population_size_; i++) {
				Chromosome<Gen>* ch = new Chromosome<Gen>(m_chromosome_size_);
				if (m_is_max_value_set_) {
					ch->set_max_value(m_max_value_);
				}
				if (m_is_min_value_set_) {
					ch->set_min_value(m_min_value_);
				}
				m_initialize_function_->run(ch);
				m_fitness_function_->run(ch);
				insert_to_population(ch);
			}
			sort_population();
			m_is_population_init_ = true;
		}

		//runs main loop
		void thread_function_() {
			size_t over_population_count, current_population_size;
			size_t children_num;
			std::queue<int> selected_parent_indexes;
			std::set<int>parents_to_execute;

			// main loop
			while (m_run_ && !m_found_min_required_solutions_ && m_generations_ < m_generation_limit_) {
				auto start = std::chrono::high_resolution_clock::now();
				std::vector<Chromosome<Gen>*> children;
				current_population_size = m_population_.size();

				//Selection phase
				auto selection_phase_start = std::chrono::high_resolution_clock::now();
				m_selection_function_->run(m_population_, selected_parent_indexes);
				size_t number_of_selected_parents = selected_parent_indexes.size();
				if (m_execute_parents_mode_ == run_t::ON) {
					for (unsigned int i = 0; i < number_of_selected_parents; i++) {
						unsigned int parent_i = selected_parent_indexes.front();
						parents_to_execute.insert(parent_i);
						selected_parent_indexes.pop();
						selected_parent_indexes.push(i);
					}
				}
				auto selection_phase_end = std::chrono::high_resolution_clock::now();
				auto selection_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(selection_phase_end - selection_phase_start);

				//Crossover phase
				auto crossover_phase_start = std::chrono::high_resolution_clock::now();
				while (!selected_parent_indexes.empty()) {
					unsigned int parent1_index = selected_parent_indexes.front();
					selected_parent_indexes.pop();
					unsigned int parent2_index = selected_parent_indexes.front();
					selected_parent_indexes.pop();
					if (prob() <= m_crossover_prob_) {
						m_crossover_function_->run(*get_kth_chromosome(parent1_index + 1), *get_kth_chromosome(parent2_index + 1), children);
					}
				}
				if (m_execute_parents_mode_ == run_t::ON) {
					while(!parents_to_execute.empty()) {
						auto it = parents_to_execute.end();
						delete *(m_population_.begin() + *(--it));
						it++;
						m_population_.erase(m_population_.begin() + *(--it));
						parents_to_execute.erase(*it);
					}
				}
				auto crossover_phase_end = std::chrono::high_resolution_clock::now();
				auto crossover_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(crossover_phase_end - crossover_phase_start);

				children_num = children.size();

				//Mutation phase
				auto mutation_phase_start = std::chrono::high_resolution_clock::now();
				for (unsigned int i = 0; i < children_num; i++) {
					if (prob() < m_chromosome_mutation_prob_) {
						m_mutation_function_->run(children.at(i));
					}
				}
				auto mutation_phase_end = std::chrono::high_resolution_clock::now();
				auto mutation_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(mutation_phase_end - mutation_phase_start);

				//Fitness phase
				auto fitness_phase_start = std::chrono::high_resolution_clock::now();
				for (unsigned int i = 0; i < children_num; i++) {
					m_fitness_function_->run(children.at(i));
				}
				auto fitness_phase_end = std::chrono::high_resolution_clock::now();
				auto fitness_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(fitness_phase_end - fitness_phase_start);

				//insertion phase
				auto insertion_phase_start = std::chrono::high_resolution_clock::now();
				size_t k = children_num;
				size_t n = m_population_.size();
				if ((k*(n + k)) / 2 <= k + (n + k)*std::log2(n + k)) { //better to insert sorted
					for (unsigned int i = 0; i < children_num; i++) {
						if (m_add_existing_chromosomes_mode_ == run_t::OFF && is_chromosome_in_population(children.at(i))) {
							delete children[i];
							continue;
						}
						insert_to_population_sorted(children.at(i));
					}
				}
				else {
					for (unsigned int i = 0; i < children_num; i++) {
						if (m_add_existing_chromosomes_mode_ == run_t::OFF && is_chromosome_in_population(children.at(i))) {
							delete children[i];
							continue;
						}
						insert_to_population(children.at(i));
					}
				}
				auto insertion_phase_end = std::chrono::high_resolution_clock::now();
				auto insertion_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(insertion_phase_end - insertion_phase_start);

				//Sort phase (if inserted not in sorted order)
				auto sort_phase_start = std::chrono::high_resolution_clock::now();
				sort_population();
				auto sort_phase_end = std::chrono::high_resolution_clock::now();
				auto sort_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(sort_phase_end - sort_phase_start);

				//Execution phase
				auto execution_phase_start = std::chrono::high_resolution_clock::now();
				if (is_over_populated(over_population_count) && m_is_execution_on_) {
					execute(over_population_count,rank_t::WEAK);
				}
				auto execution_phase_end = std::chrono::high_resolution_clock::now();
				auto execution_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(execution_phase_end - execution_phase_start);

				m_generations_++;

				//Solution Phase
				auto solution_phase_start = std::chrono::high_resolution_clock::now();
				while (m_number_of_found_solutions_ < m_min_number_of_solutions_ && m_population_.at(m_number_of_found_solutions_)->get_fitness() >= m_min_solution_fitness_) {
					if (!solution_exist(m_population_.at(m_number_of_found_solutions_))) { // solution doesn't exist
						m_solutions_.push_back(m_population_.at(m_number_of_found_solutions_));
						m_number_of_found_solutions_++;
					}
					else { //solution exists -> remove it
						remove_by_index(m_number_of_found_solutions_);
					}
				}
				auto solution_phase_end = std::chrono::high_resolution_clock::now();
				auto solution_phase_duration = std::chrono::duration_cast<std::chrono::microseconds>(solution_phase_end - solution_phase_start);

				if (m_number_of_found_solutions_ >= m_min_number_of_solutions_) {
					m_found_min_required_solutions_ = true;
				}
				current_population_size = m_population_.size();

				//Fitness mean phase
				double mean, sum = 0;
				for (unsigned int i = 0; i < current_population_size; i++) {
					sum += m_population_.at(i)->get_fitness();
				}
				mean = sum / current_population_size;
				m_mean_diff_ = mean - m_mean_;
				m_mean_ = mean;
				//Fitness median,highest,lowest phase
				double median = m_population_.at((unsigned int)std::ceil(current_population_size / 2))->get_fitness();
				double highest_fitness = m_population_[0]->get_fitness();
				double lowest_fitness = m_population_[m_population_.size() - 1]->get_fitness();
				auto stop = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

				//Print phase
				if (m_generations_ - m_previous_print_generation_ >= m_print_after_passed_generations_)
				{
					m_previous_print_generation_ = m_generations_;

					std::cout.precision(10);
					std::cout << "Summary Iteration [" << m_generations_ << "]\n Selection time: \t" << selection_phase_duration.count() << "\tmicroseconds\n Crossover time: \t" << crossover_phase_duration.count() <<
						"\tmicroseconds\n Mutation time: \t" << mutation_phase_duration.count() << "\tmicroseconds\n Fitness time: \t\t" << fitness_phase_duration.count() << "\tmicroseconds\n Insertion time: \t" << insertion_phase_duration.count() <<
						"\tmicroseconds\n Sort time: \t\t" << sort_phase_duration.count() << "\tmicroseconds\n Execution time: \t" << execution_phase_duration.count() << "\tmicroseconds\n Solution time: \t" << solution_phase_duration.count() << "\tmicroseconds\n Total time: \t\t" << duration.count() << "\tmicroseconds\n Current population size: \t" <<
						m_population_.size() << "\n New chromosomes: \t\t" << children_num << "\n Executed chromosomes: \t\t" << over_population_count <<"\n Top fitness: \t\t\t"<< highest_fitness << "\n Low fitness: \t\t\t" << lowest_fitness <<"\n Fitness Mean: \t\t\t" << mean << "\n Mean diff: \t\t\t" << m_mean_diff_;
					std::cout << "\n Median Fitness: \t\t" << median <<"\n Found solutions: \t\t[" << m_number_of_found_solutions_ << "/" << m_min_number_of_solutions_ << "]\n-------------------------------------------------\n";
				}
				//Save mean phase
				if (m_generations_ - m_previous_saved_mean_generation_ >= m_save_mean_after_passed_generations_)
				{
					m_previous_saved_mean_generation_ = m_generations_;
					m_saved_generation_.push_back(m_generations_);
					m_saved_highest_fitness_.push_back(highest_fitness);
					m_saved_lowest_fitness_.push_back(lowest_fitness);
					m_saved_fitness_mean_.push_back(m_mean_);
				}
				//Eradicate if required
				if (should_eradicate() && prob() <= m_eradicate_prob_) {
					eradicate(m_eradicate_percentage_);
					init_population_();
					m_num_of_eradications_++;
				}
			} //while true
			if (m_found_min_required_solutions_) {
				std::cout << "Done !!!\n";
			}
		}

		//returns if population should be eradicated
		bool should_eradicate() {
			if (m_generations_ == 100) {
				printf("");
			}
			if (m_previous_recorded_mean_ == NULL) {
				m_previous_recorded_mean_ = m_mean_;
				return false;
			}
			if (m_eradicate_generations_ <= 0) {
				return false;
			}
			if (m_generations_ - m_previous_recorded_generation_ >= m_eradicate_generations_) {
				// fitness mean should only go up by definition, so slope is required to be > 0
				double progress_slope = (m_mean_ - m_previous_recorded_mean_) / (m_generations_ - m_previous_recorded_generation_);
				m_previous_recorded_generation_ = m_generations_;
				m_previous_recorded_mean_ = m_mean_;
				return progress_slope < m_eradicate_slope_;
			}
			return false;
		}
	};

	namespace selection_func_ns
	{
		template <typename Gen>
		class ExploreExploitSelection : public SelectionFunction<Gen> {
		public:
			void run(const std::vector<Chromosome<Gen>*>& population, std::queue<int>& population_index_selected) override {
				size_t population_size = population.size();
				double number_of_selections = std::ceil(population_size * SelectionFunction<Gen>::pm_select_ratio) - 1;
				unsigned int selections_made = 0;

				//Clean queue if recieved full
				while (!population_index_selected.empty()) {
					population_index_selected.pop();
				}

				//Crossover with himself
				if (population_size == 1) {
					population_index_selected.push(0);
					population_index_selected.push(0);
					return;
				}

				//we reach here if population size is atleast 2, max number of selections is population size - 1
				for (unsigned int i = 0; i < population_size; i++) {

					if (selections_made >= number_of_selections) {
						break;
					}
				//exploit
					population_index_selected.push(i);
					if (prob() <= SelectionFunction<Gen>::pm_exploit_prob) {
						population_index_selected.push(i + 1);
					}
				//explore
					else {
						population_index_selected.push(random(0, population_size-1));
					}
					selections_made++;
				}
			}
		};
	}

	namespace crossover_func_ns
	{
		template<typename Gen>
		class TwoPointCrossover : public CrossoverFunction<Gen> {
		private:
			unsigned int m_index_;
		public:
			TwoPointCrossover(unsigned int crossover_index) {
				m_index_ = crossover_index;
			}
			unsigned int get_crossover_index() {
				return m_index_;
			}
			void set_crossover_index(unsigned int crossover_index) {
				m_index_ = crossover_index;
			}
			void run(const Chromosome<Gen>& chromsome1, const Chromosome<Gen>& chromsome2, std::vector<Chromosome<Gen>*>& children) const override {
				if (prob() > CrossoverFunction<Gen>::pm_crossover_prob) {
					return;
				}
				size_t size = chromsome1.get_size();
				if (size != chromsome2.get_size()) {
					return;
				}
				if (chromsome1.get_max_value() != chromsome2.get_max_value() || chromsome1.get_min_value() != chromsome2.get_min_value()) {
					return;
				}
				Chromosome<Gen>* new_ch1 = new Chromosome<Gen>(size);
				Chromosome<Gen>* new_ch2 = new Chromosome<Gen>(size);
				if (m_index_ >= size) {
					*new_ch1 = chromsome1;
					*new_ch2 = chromsome2;
					children.push_back(new_ch1);
					children.push_back(new_ch2);
					return;
				}
				new_ch1->set_max_value(chromsome1.get_max_value());
				new_ch1->set_min_value(chromsome1.get_min_value());
				new_ch2->set_max_value(chromsome1.get_max_value());
				new_ch2->set_min_value(chromsome1.get_min_value());
				for (unsigned int i = 0; i < m_index_; i++) {
					new_ch1->init_push(chromsome1.get_gene_at(i));
					new_ch2->init_push(chromsome2.get_gene_at(i));
				}
				for (unsigned int i = m_index_; i < size; i++) {
					new_ch1->init_push(chromsome2.get_gene_at(i));
					new_ch2->init_push(chromsome1.get_gene_at(i));
				}
				new_ch1->set_is_init();
				new_ch2->set_is_init();
				children.push_back(new_ch1);
				children.push_back(new_ch2);
			}
		};
	}
}
#endif // !GEN_LIB_H