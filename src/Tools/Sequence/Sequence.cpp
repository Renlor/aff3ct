#include <set>
#include <thread>
#include <utility>
#include <sstream>
#include <fstream>
#include <cstring>
#include <numeric>
#include <exception>
#include <algorithm>

#include "Tools/Display/rang_format/rang_format.h"
#include "Tools/Thread_pinning/Thread_pinning.hpp"
#include "Tools/Exception/exception.hpp"
#include "Module/Module.hpp"
#include "Module/Task.hpp"
#include "Module/Socket.hpp"
#include "Module/Probe/Probe.hpp"
#include "Module/Loop/Loop.hpp"
#include "Module/Router/Router.hpp"
#include "Module/Adaptor/Adaptor.hpp"
#include "Tools/Sequence/Sequence.hpp"

using namespace aff3ct;
using namespace aff3ct::tools;

Sequence
::Sequence(const std::vector<const module::Task*> &firsts,
           const std::vector<const module::Task*> &lasts,
           const std::vector<const module::Task*> &exclusions,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids)
: n_threads(n_threads),
  sequences(n_threads, nullptr),
  modules(n_threads),
  all_modules(n_threads),
  mtx_exception(new std::mutex()),
  force_exit_loop(new std::atomic<bool>(false)),
  tasks_inplace(false),
  thread_pinning(thread_pinning),
  puids(puids),
  no_copy_mode(true),
  saved_exclusions(exclusions)
{
#ifndef AFF3CT_HWLOC
	if (thread_pinning)
		std::clog << rang::tag::warning << "AFF3CT has not been linked with the 'hwloc' library, the 'thread_pinning' "
		                                   "option of the 'tools::Sequence' will have no effect." << std::endl;
#endif

	if (thread_pinning && puids.size() < n_threads)
	{
		std::stringstream message;
		message << "'puids.size()' has to be greater or equal to 'n_threads' ('puids.size()' = " << puids.size()
		        << " , 'n_threads' = " << n_threads << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	this->init<tools::Sub_sequence_const,const module::Task>(firsts, lasts, exclusions);
}

Sequence
::Sequence(const std::vector<const module::Task*> &firsts,
           const std::vector<const module::Task*> &lasts,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids)
: Sequence(firsts, lasts, {}, n_threads, thread_pinning, puids)
{
}

Sequence
::Sequence(const std::vector<const module::Task*> &firsts,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids)
: Sequence(firsts, {}, {}, n_threads, thread_pinning, puids)
{
}

Sequence
::Sequence(const module::Task &first,
           const module::Task &last,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids)
: Sequence({&first}, {&last}, n_threads, thread_pinning, puids)
{
}

Sequence
::Sequence(const module::Task &first,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids)
: Sequence({&first}, n_threads, thread_pinning, puids)
{
}

std::vector<const module::Task*> exclusions_convert_to_const(const std::vector<module::Task*> &exclusions)
{
	std::vector<const module::Task*> exclusions_const;
	for (auto exception : exclusions)
		exclusions_const.push_back(exception);
	return exclusions_const;
}

Sequence
::Sequence(const std::vector<module::Task*> &firsts,
           const std::vector<module::Task*> &lasts,
           const std::vector<module::Task*> &exclusions,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids,
           const bool tasks_inplace)
: n_threads(n_threads),
  sequences(n_threads, nullptr),
  modules(tasks_inplace ? n_threads -1 : n_threads),
  all_modules(n_threads),
  mtx_exception(new std::mutex()),
  force_exit_loop(new std::atomic<bool>(false)),
  tasks_inplace(tasks_inplace),
  thread_pinning(thread_pinning),
  puids(puids),
  no_copy_mode(true),
  saved_exclusions(exclusions_convert_to_const(exclusions))
{
	if (thread_pinning && puids.size() < n_threads)
	{
		std::stringstream message;
		message << "'puids.size()' has greater or equal to 'n_threads' ('puids.size()' = " << puids.size()
		        << " , 'n_threads' = " << n_threads << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	if (tasks_inplace)
		this->init<tools::Sub_sequence,module::Task>(firsts, lasts, exclusions);
	else
	{
		std::vector<const module::Task*> firsts_bis;
		for (auto first : firsts) firsts_bis.push_back(first);
		std::vector<const module::Task*> lasts_bis;
		for (auto last : lasts) lasts_bis.push_back(last);
		std::vector<const module::Task*> exclusions_bis;
		for (auto exception : exclusions) exclusions_bis.push_back(exception);
		this->init<tools::Sub_sequence_const,const module::Task>(firsts_bis, lasts_bis, exclusions_bis);
	}
}

Sequence
::Sequence(const std::vector<module::Task*> &firsts,
           const std::vector<module::Task*> &lasts,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids,
           const bool tasks_inplace)
: Sequence(firsts, lasts, {}, n_threads, thread_pinning, puids, tasks_inplace)
{
}

Sequence
::Sequence(const std::vector<module::Task*> &firsts,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids,
           const bool tasks_inplace)
: Sequence(firsts, {}, {}, n_threads, thread_pinning, puids, tasks_inplace)
{
}

Sequence
::Sequence(module::Task &first,
           module::Task &last,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids,
           const bool tasks_inplace)
: Sequence({&first}, {&last}, n_threads, thread_pinning, puids, tasks_inplace)
{
}

Sequence
::Sequence(module::Task &first,
           const size_t n_threads,
           const bool thread_pinning,
           const std::vector<size_t> &puids,
           const bool tasks_inplace)
: Sequence({&first}, n_threads, thread_pinning, puids, tasks_inplace)
{
}

Sequence
::~Sequence()
{
	for (auto s : this->sequences)
		this->delete_tree(s);
}

template <class SS, class TA>
void Sequence
::init(const std::vector<TA*> &firsts, const std::vector<TA*> &lasts, const std::vector<TA*> &exclusions)
{
	if (this->is_thread_pinning())
		Thread_pinning::pin(this->puids[0]);

	if (firsts.size() == 0)
	{
		std::stringstream message;
		message << "'firsts.size()' has to be strictly greater than 0.";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	if (this->n_threads == 0)
	{
		std::stringstream message;
		message << "'n_threads' has to be strictly greater than 0.";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	for (auto exclusion : exclusions)
	{
		if (std::find(firsts.begin(), firsts.end(), exclusion) != firsts.end())
		{
			std::stringstream message;
			message << "'exclusion' can't be contained in the 'firsts' vector ("
			        << "'exclusion'"                 << " = " << +exclusion            << ", "
			        << "'exclusion->get_name()'"     << " = " << exclusion->get_name() << ").";
			throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
		}

		if (std::find(lasts.begin(), lasts.end(), exclusion) != lasts.end())
		{
			std::stringstream message;
			message << "'exclusion' can't be contained in the 'lasts' vector ("
			        << "'exclusion'"                 << " = " << +exclusion            << ", "
			        << "'exclusion->get_name()'"     << " = " << exclusion->get_name() << ").";
			throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
		}
	}

	auto root = new Generic_node<SS>(nullptr, {}, nullptr, 0, 0, 0);
	root->set_contents(nullptr);
	size_t ssid = 0, taid = 0;
	std::vector<TA*> loops;
	std::vector<TA*> real_lasts;

	this->lasts_tasks_id.clear();
	this->firsts_tasks_id.clear();
	auto last_subseq = root;
	for (auto first : firsts)
	{
		auto contents = last_subseq->get_contents();
		this->firsts_tasks_id.push_back(contents ? contents->tasks_id[contents->tasks_id.size() -1] : 0);
		last_subseq = this->init_recursive<SS,TA>(last_subseq,
		                                          ssid,
		                                          taid,
		                                          loops,
		                                          *first,
		                                          *first,
		                                          lasts,
		                                          exclusions,
		                                          this->lasts_tasks_id,
		                                          real_lasts);
	}

	std::stringstream real_lasts_ss;
	for (size_t rl = 0; rl < real_lasts.size(); rl++)
		real_lasts_ss << "'real_lasts" << "[" << rl << "]'"             << " = " << +real_lasts[rl] << ", "
		              << "'real_lasts" << "[" << rl << "]->get_name()'" << " = " <<  real_lasts[rl]->get_name()
		              << ((rl < real_lasts.size() -1) ? ", " : "");

	for (auto last : lasts)
	{
		if (std::find(real_lasts.begin(), real_lasts.end(), last) == real_lasts.end())
		{
			std::stringstream message;
			message << "'last' is not contained in the 'real_lasts[" << real_lasts.size() << "]' vector ("
			        << "'last'"                 << " = " << +last            << ", "
			        << "'last->get_name()'"     << " = " << last->get_name() << ", "
			        << real_lasts_ss.str() << ").";
			throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
		}
	}

	this->n_tasks = taid;
	this->_init<SS>(root);
	this->update_firsts_and_lasts_tasks();
	this->gen_processes();
}

Sequence* Sequence
::clone() const
{
	auto c = new Sequence(*this);

	c->tasks_inplace = false;
	c->modules.resize(c->get_n_threads());

	std::vector<const module::Task*> firsts_tasks;
	for (auto ta : this->get_firsts_tasks()[0])
		firsts_tasks.push_back(ta);

	std::vector<const module::Task*> lasts_tasks;
	for (auto ta : this->get_lasts_tasks()[0])
		lasts_tasks.push_back(ta);

	c->init<tools::Sub_sequence_const,const module::Task>(firsts_tasks, lasts_tasks, this->saved_exclusions);
	c->mtx_exception.reset(new std::mutex());
	c->force_exit_loop.reset(new std::atomic<bool>(false));
	return c;
}

void Sequence
::set_thread_pinning(const bool thread_pinning, const std::vector<size_t> &puids)
{
	if (thread_pinning && puids.size() < n_threads)
	{
		std::stringstream message;
		message << "'puids.size()' has greater or equal to 'n_threads' ('puids.size()' = " << puids.size()
		        << " , 'n_threads' = " << n_threads << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	this->thread_pinning = thread_pinning;
	this->puids = puids;
}

bool Sequence
::is_thread_pinning()
{
	return this->thread_pinning;
}

std::vector<std::vector<module::Module*>> Sequence
::get_modules_per_threads() const
{
	std::vector<std::vector<module::Module*>> modules_per_threads(this->all_modules.size());
	size_t tid = 0;
	for (auto &e : this->all_modules)
	{
		for (auto &ee : e)
			modules_per_threads[tid].push_back(ee);
		tid++;
	}
	return modules_per_threads;
}

std::vector<std::vector<module::Module*>> Sequence
::get_modules_per_types() const
{
	std::vector<std::vector<module::Module*>> modules_per_types(this->all_modules[0].size());
	for (auto &e : this->all_modules)
	{
		size_t mid = 0;
		for (auto &ee : e)
			modules_per_types[mid++].push_back(ee);
	}
	return modules_per_types;
}

std::vector<std::vector<module::Task*>> Sequence
::get_tasks_per_threads() const
{
	std::vector<std::vector<module::Task*>> tasks_per_threads(this->n_threads);

	std::function<void(Generic_node<Sub_sequence>*, const size_t)> get_tasks_recursive =
		[&](Generic_node<Sub_sequence>* cur_ss, const size_t tid)
		{
			tasks_per_threads[tid].insert(tasks_per_threads[tid].end(),
			                              cur_ss->get_c()->tasks.begin(),
			                              cur_ss->get_c()->tasks.end());

			for (auto c : cur_ss->get_children())
				get_tasks_recursive(c, tid);
		};

	for (size_t tid = 0; tid < this->n_threads; tid++)
		get_tasks_recursive(this->sequences[tid], tid);

	return tasks_per_threads;
}

std::vector<std::vector<module::Task*>> Sequence
::get_tasks_per_types() const
{
	std::vector<std::vector<module::Task*>> tasks_per_types(this->n_tasks);

	std::function<void(Generic_node<Sub_sequence>*, size_t&)> get_tasks_recursive =
		[&](Generic_node<Sub_sequence>* cur_ss, size_t &mid)
		{
			for (auto &t : cur_ss->get_c()->tasks)
				tasks_per_types[mid++].push_back(t);

			for (auto c : cur_ss->get_children())
				get_tasks_recursive(c, mid);
		};

	for (size_t tid = 0; tid < this->n_threads; tid++)
	{
		size_t mid = 0;
		get_tasks_recursive(this->sequences[tid], mid);
	}

	return tasks_per_types;
}

void Sequence
::_exec(const size_t tid,
        std::function<bool(const std::vector<int>&)> &stop_condition,
        Generic_node<Sub_sequence>* sequence)
{
	if (this->is_thread_pinning())
		Thread_pinning::pin(this->puids[tid]);

	std::function<void(Generic_node<Sub_sequence>*, std::vector<int>&)> exec_sequence =
		[&exec_sequence](Generic_node<Sub_sequence>* cur_ss, std::vector<int>& statuses)
		{
			auto type = cur_ss->get_c()->type;
			auto &tasks = cur_ss->get_c()->tasks;
			auto &tasks_id = cur_ss->get_c()->tasks_id;
			auto &processes = cur_ss->get_c()->processes;

			if (type == subseq_t::LOOP)
			{
				while (!(statuses[tasks_id[0]] = processes[0]()))
					exec_sequence(cur_ss->get_children()[0], statuses);
				static_cast<module::Loop&>(tasks[0]->get_module()).reset();
				exec_sequence(cur_ss->get_children()[1], statuses);
			}
			else
			{
				for (size_t p = 0; p < processes.size(); p++)
					statuses[tasks_id[p]] = processes[p]();
				for (auto c : cur_ss->get_children())
					exec_sequence(c, statuses);
			}
		};

	std::vector<int> statuses(this->n_tasks, 0);
	try
	{
		do
		{
			std::fill(statuses.begin(), statuses.end(), module::status_t::SKIPPED);
			try
			{
				exec_sequence(sequence, statuses);
			}
			catch (tools::processing_aborted const&)
			{
				// do nothing, this is normal
			}
		}
		while (!*force_exit_loop && !stop_condition(statuses));
	}
	catch (tools::waiting_canceled const&)
	{
		// do nothing, this is normal
	}
	catch (std::exception const& e)
	{
		*force_exit_loop = true;

		this->mtx_exception->lock();

		auto save = tools::exception::no_backtrace;
		tools::exception::no_backtrace = true;
		std::string msg = e.what(); // get only the function signature
		tools::exception::no_backtrace = save;

		if (std::find(this->prev_exception_messages.begin(), this->prev_exception_messages.end(), msg) ==
			this->prev_exception_messages.end())
		{
			this->prev_exception_messages.push_back(msg); // save only the function signature
			this->prev_exception_messages_to_display.push_back(e.what()); // with backtrace if debug mode
		}
		else if (std::strlen(e.what()) > this->prev_exception_messages_to_display.back().size())
			this->prev_exception_messages_to_display[prev_exception_messages_to_display.size() -1] = e.what();

		this->mtx_exception->unlock();
	}

	if (this->is_thread_pinning())
		Thread_pinning::unpin();
}

void Sequence
::_exec_without_statuses(const size_t tid,
                         std::function<bool()> &stop_condition,
                         Generic_node<Sub_sequence>* sequence)
{
	if (this->is_thread_pinning())
		Thread_pinning::pin(this->puids[tid]);

	const bool no_copy_mode = this->is_no_copy_mode();
	std::function<void(Generic_node<Sub_sequence>*)> exec_sequence =
		[&exec_sequence, no_copy_mode](Generic_node<Sub_sequence>* cur_ss)
		{
			auto type = cur_ss->get_c()->type;
			auto &tasks = cur_ss->get_c()->tasks;
			auto &processes = cur_ss->get_c()->processes;

			if (type == subseq_t::LOOP)
			{
				while (!processes[0]())
					exec_sequence(cur_ss->get_children()[0]);
				static_cast<module::Loop&>(tasks[0]->get_module()).reset();
				exec_sequence(cur_ss->get_children()[1]);
			}
			else
			{
				for (auto &process : processes)
					process();
				for (auto c : cur_ss->get_children())
					exec_sequence(c);
			}
		};

	try
	{
		do
		{
			try
			{
				exec_sequence(sequence);
			}
			catch (tools::processing_aborted const&)
			{
				// do nothing, this is normal
			}
		}
		while (!*force_exit_loop && !stop_condition());
	}
	catch (tools::waiting_canceled const&)
	{
		// do nothing, this is normal
	}
	catch (std::exception const& e)
	{
		*force_exit_loop = true;

		this->mtx_exception->lock();

		auto save = tools::exception::no_backtrace;
		tools::exception::no_backtrace = true;
		std::string msg = e.what(); // get only the function signature
		tools::exception::no_backtrace = save;

		if (std::find(this->prev_exception_messages.begin(), this->prev_exception_messages.end(), msg) ==
			this->prev_exception_messages.end())
		{
			this->prev_exception_messages.push_back(msg); // save only the function signature
			this->prev_exception_messages_to_display.push_back(e.what()); // with backtrace if debug mode
		}
		else if (std::strlen(e.what()) > this->prev_exception_messages_to_display.back().size())
			this->prev_exception_messages_to_display[prev_exception_messages_to_display.size() -1] = e.what();

		this->mtx_exception->unlock();
	}

	if (this->is_thread_pinning())
		Thread_pinning::unpin();
}

void Sequence
::exec(std::function<bool(const std::vector<int>&)> stop_condition)
{
	if (this->is_no_copy_mode())
		this->gen_processes(true);

	std::vector<std::thread> threads(n_threads);
	for (size_t tid = 1; tid < n_threads; tid++)
		threads[tid] = std::thread(&Sequence::_exec, this, tid, std::ref(stop_condition), std::ref(this->sequences[tid]));

	this->_exec(0, stop_condition, this->sequences[0]);

	for (size_t tid = 1; tid < n_threads; tid++)
		threads[tid].join();

	if (this->is_no_copy_mode())
	{
		this->reset_no_copy_mode();
		this->gen_processes(false);
	}

	if (!this->prev_exception_messages_to_display.empty())
	{
		*force_exit_loop = false;
		throw std::runtime_error(this->prev_exception_messages_to_display.back());
	}
}

void Sequence
::exec(std::function<bool()> stop_condition)
{
	if (this->is_no_copy_mode())
		this->gen_processes(true);

	std::vector<std::thread> threads(n_threads);
	for (size_t tid = 1; tid < n_threads; tid++)
	{
		threads[tid] = std::thread(&Sequence::_exec_without_statuses, this, tid, std::ref(stop_condition),
		                           std::ref(this->sequences[tid]));
	}

	this->_exec_without_statuses(0, stop_condition, this->sequences[0]);

	for (size_t tid = 1; tid < n_threads; tid++)
		threads[tid].join();

	if (this->is_no_copy_mode())
	{
		this->reset_no_copy_mode();
		this->gen_processes(false);
	}

	if (!this->prev_exception_messages_to_display.empty())
	{
		*force_exit_loop = false;
		throw std::runtime_error(this->prev_exception_messages_to_display.back());
	}
}

int Sequence
::exec(const size_t tid)
{
	if (tid >= this->sequences.size())
	{
		std::stringstream message;
		message << "'tid' has to be smaller than 'sequences.size()' ('tid' = " << tid
		        << ", 'sequences.size()' = " << this->sequences.size() << ").";
		throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
	}

	std::function<void(Generic_node<Sub_sequence>*, int&)> exec_sequence =
		[&exec_sequence](Generic_node<Sub_sequence>* cur_ss, int& ret)
		{
			auto type = cur_ss->get_c()->type;
			auto &tasks = cur_ss->get_c()->tasks;
			if (type == subseq_t::LOOP)
			{
				while (!tasks[0]->exec())
					exec_sequence(cur_ss->get_children()[0], ret);
				ret++;
				static_cast<module::Loop&>(tasks[0]->get_module()).reset();
				exec_sequence(cur_ss->get_children()[1], ret);
			}
			else
			{
				for (size_t ta = 0; ta < tasks.size(); ta++)
					ret += tasks[ta]->exec();
				for (auto c : cur_ss->get_children())
					exec_sequence(c, ret);
			}
		};

	int ret = 0;
	exec_sequence(this->sequences[tid], ret);
	return ret;
}

template <class SS, class TA>
Generic_node<SS>* Sequence
::init_recursive(Generic_node<SS> *cur_subseq,
                 size_t &ssid,
                 size_t &taid,
                 std::vector<TA*> &loops,
                 TA &first,
                 TA &current_task,
                 const std::vector<TA*> &lasts,
                 const std::vector<TA*> &exclusions,
                 std::vector<size_t> &real_lasts_id,
                 std::vector<TA*> &real_lasts)
{
	if (this->tasks_inplace && !current_task.is_autoalloc())
	{
		std::stringstream message;
		message << "When 'tasks_inplace' is set to true 'current_task' should be in autoalloc mode ("
		        << "'current_task.get_name()'"              << " = " << current_task.get_name()              << ", "
		        << "'current_task.get_module().get_name()'" << " = " << current_task.get_module().get_name() << ").";
		throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
	}

	if (dynamic_cast<const module::Adaptor*>(&current_task.get_module()) && !this->tasks_inplace)
	{
		std::stringstream message;
		message << "'module::Adaptor' objects are not supported if 'tasks_inplace' is set to false.";
		throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
	}

	auto it = std::find(real_lasts.begin(), real_lasts.end(), &current_task);
	if (it != real_lasts.end())
	{
		real_lasts.erase(it);
		auto dist = std::distance(real_lasts.begin(), it);
		real_lasts_id.erase(real_lasts_id.begin() + dist);
	}

	if (auto loop = dynamic_cast<const module::Loop*>(&current_task.get_module()))
	{
		if (std::find(loops.begin(), loops.end(), &current_task) == loops.end())
		{
			loops.push_back(&current_task);
			Generic_node<SS>* node_loop = nullptr;
			if (&first == &current_task)
				node_loop = cur_subseq;
			else
			{
				ssid++;
				node_loop = new Generic_node<SS>(cur_subseq, {}, nullptr, cur_subseq->get_depth() +1, 0, 0);
			}

			auto node_loop_son0 = new Generic_node<SS>(node_loop, {}, nullptr, node_loop->get_depth() +1, 0, 0);
			auto node_loop_son1 = new Generic_node<SS>(node_loop, {}, nullptr, node_loop->get_depth() +1, 0, 1);
			node_loop->add_child(node_loop_son0);
			node_loop->add_child(node_loop_son1);

			node_loop->set_contents(new SS());
			node_loop_son0->set_contents(new SS());
			node_loop_son1->set_contents(new SS());

			node_loop->get_c()->tasks.push_back(&current_task);
			node_loop->get_c()->tasks_id.push_back(taid++);
			node_loop->get_c()->type = subseq_t::LOOP;
			node_loop->get_c()->id = ssid++;

			if (!cur_subseq->get_children().size())
				cur_subseq->add_child(node_loop);

			if (loop->tasks[0]->sockets[2]->get_bound_sockets().size() == 1)
			{
				node_loop_son0->get_c()->id = ssid++;
				auto &t = loop->tasks[0]->sockets[2]->get_bound_sockets()[0]->get_task();
				Sequence::init_recursive<SS,TA>(node_loop_son0, ssid, taid, loops, first, t, lasts, exclusions,
				                                real_lasts_id, real_lasts);
			}
			else
			{
				std::stringstream message;
				message << "'loop->tasks[0]->sockets[2]->get_bound_sockets().size()' has to be equal to 1 ("
				        << "'loop->tasks[0]->sockets[2]->get_bound_sockets().size()' = "
				        << loop->tasks[0]->sockets[2]->get_bound_sockets().size() << ").";
				throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
			}

			if (loop->tasks[0]->sockets[3]->get_bound_sockets().size() == 1)
			{
				node_loop_son1->get_c()->id = ssid++;
				auto &t = loop->tasks[0]->sockets[3]->get_bound_sockets()[0]->get_task();
				return Sequence::init_recursive<SS,TA>(node_loop_son1, ssid, taid, loops, first, t, lasts, exclusions,
				                                       real_lasts_id, real_lasts);
			}
			else
			{
				std::stringstream message;
				message << "'loop->tasks[0]->sockets[3]->get_bound_sockets().size()' has to be equal to 1 ("
				        << "'loop->tasks[0]->sockets[3]->get_bound_sockets().size()' = "
				        << loop->tasks[0]->sockets[3]->get_bound_sockets().size() << ").";
				throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
			}
		}
	}
	else
	{
		if (cur_subseq->get_contents() == nullptr)
			cur_subseq->set_contents(new SS());

		cur_subseq->get_c()->tasks.push_back(&current_task);
		cur_subseq->get_c()->tasks_id.push_back(taid++);

		bool is_last = true;
		Generic_node<SS>* last_subseq = nullptr;
		if (std::find(lasts.begin(), lasts.end(), &current_task) == lasts.end())
		{
			for (auto &s : current_task.sockets)
			{
				if (current_task.get_socket_type(*s) == module::socket_t::SIN_SOUT ||
					current_task.get_socket_type(*s) == module::socket_t::SOUT)
				{
					auto bss = s->get_bound_sockets();
					for (auto &bs : bss)
					{
						if (bs != nullptr)
						{
							auto &t = bs->get_task();
							if (std::find(exclusions.begin(), exclusions.end(), &t) == exclusions.end())
							{
								if (t.is_last_input_socket(*bs) || dynamic_cast<const module::Loop*>(&t.get_module()))
								{
									is_last = false;
									last_subseq = Sequence::init_recursive<SS,TA>(cur_subseq, ssid, taid, loops, first,
									                                              t, lasts, exclusions, real_lasts_id,
									                                              real_lasts);
								}
							}
						}
					}
				}
				else if (current_task.get_socket_type(*s) == module::socket_t::SIN)
				{
					if (s->get_bound_sockets().size() > 1)
					{
						std::stringstream message;
						message << "'s->get_bound_sockets().size()' has to be smaller or equal to 1 ("
						        << "'s->get_bound_sockets().size()'"         << " = " << s->get_bound_sockets().size() << ", "
						        << "'get_socket_type(*s)'"                   << " = " << "socket_t::SIN"               << ", "
						        << "'s->get_name()'"                         << " = " << s->get_name()                 << ", "
						        << "'s->get_task().get_name()'"              << " = " << s->get_task().get_name()      << ", "
						        << "'s->get_task().get_module().get_name()'" << " = " << s->get_task().get_module().get_name()
						        << ").";
						throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
					}
				}
			}
		}

		if (is_last && std::find(real_lasts.begin(), real_lasts.end(), &current_task) == real_lasts.end())
		{
			real_lasts.push_back(&current_task);
			real_lasts_id.push_back(cur_subseq->get_contents()->tasks_id.back());
		}
		if (last_subseq)
			return last_subseq;
	}

	return cur_subseq;
}

template <class SS, class MO>
void Sequence
::duplicate(const Generic_node<SS> *sequence)
{
	std::set<MO*> modules_set;

	std::function<void(const Generic_node<SS>*)> collect_modules_list;
	collect_modules_list = [&](const Generic_node<SS> *node)
	{
		if (node != nullptr)
		{
			if (node->get_c())
				for (auto ta : node->get_c()->tasks)
					modules_set.insert(&ta->get_module());
			for (auto c : node->get_children())
				collect_modules_list(c);
		}
	};
	collect_modules_list(sequence);

	std::vector<MO*> modules_vec;
	for (auto m : modules_set)
		modules_vec.push_back(m);

	// clone the modules
	for (size_t tid = 0; tid < this->n_threads - (this->tasks_inplace ? 1 : 0); tid++)
	{
		if (this->is_thread_pinning())
		{
			const auto real_tid = tid + (this->tasks_inplace ? 1 : 0);
			Thread_pinning::pin(this->puids[real_tid]);
		}

		this->modules[tid].resize(modules_vec.size());
		this->all_modules[tid + (this->tasks_inplace ? 1 : 0)].resize(modules_vec.size());
		for (size_t m = 0; m < modules_vec.size(); m++)
		{
			try
			{
				this->modules[tid][m].reset(modules_vec[m]->clone());
			}
			catch (std::exception &e)
			{
				std::cerr << rang::tag::error << "Module clone failed when trying to duplicate the sequence: module "
				                              << "name is '" << modules_vec[m]->get_name() << "'." << std::endl;
				throw e;
			}
			this->all_modules[tid + (this->tasks_inplace ? 1 : 0)][m] = this->modules[tid][m].get();
		}

		if (this->is_thread_pinning())
			Thread_pinning::unpin();
	}

	auto get_module_id = [](const std::vector<MO*> &modules, const module::Module &module)
	{
		int m_id;
		for (m_id = 0; m_id < (int)modules.size(); m_id++)
			if (modules[m_id] == &module)
				return m_id;
		return -1;
	};

	auto get_task_id = [](const std::vector<std::shared_ptr<module::Task>> &tasks, const module::Task &task)
	{
		int t_id;
		for (t_id = 0; t_id < (int)tasks.size(); t_id++)
			if (tasks[t_id].get() == &task)
				return t_id;
		return -1;
	};

	auto get_socket_id = [](const std::vector<std::shared_ptr<module::Socket>> &sockets, const module::Socket &socket)
	{
		int s_id;
		for (s_id = 0; s_id < (int)sockets.size(); s_id++)
			if (sockets[s_id].get() == &socket)
				return s_id;
		return -1;
	};

	std::function<void(const Generic_node<SS>*,
	                         Generic_node<Sub_sequence>*,
	                   const size_t)> duplicate_sequence;

	duplicate_sequence = [&](const Generic_node<SS>           *sequence_ref,
	                               Generic_node<Sub_sequence> *sequence_cpy,
	                         const size_t thread_id)
	{
		if (sequence_ref != nullptr && sequence_ref->get_c())
		{
			auto ss_ref = sequence_ref->get_c();
			auto ss_cpy = new Sub_sequence();

			ss_cpy->type = ss_ref->type;
			ss_cpy->id = ss_ref->id;
			ss_cpy->tasks_id = ss_ref->tasks_id;
			for (auto t_ref : ss_ref->tasks)
			{
				auto &m_ref = t_ref->get_module();

				auto m_id = get_module_id(modules_vec, m_ref);
				auto t_id = get_task_id(m_ref.tasks, *t_ref);

				assert(m_id != -1);
				assert(t_id != -1);

				// add the task to the sub-sequence
				ss_cpy->tasks.push_back(this->all_modules[thread_id][m_id]->tasks[t_id].get());

				// replicate the sockets binding
				for (size_t s_id = 0; s_id < t_ref->sockets.size(); s_id++)
				{
					if (t_ref->get_socket_type(*t_ref->sockets[s_id]) == module::socket_t::SIN_SOUT ||
					    t_ref->get_socket_type(*t_ref->sockets[s_id]) == module::socket_t::SIN)
					{
						const module::Socket* s_ref_out = nullptr;
						try { s_ref_out = &t_ref->sockets[s_id]->get_bound_socket(); } catch (...) {}
						if (s_ref_out)
						{
							auto &t_ref_out = s_ref_out->get_task();
							auto &m_ref_out = t_ref_out.get_module();

							auto m_id_out = get_module_id(modules_vec, m_ref_out);

							if (m_id_out != -1)
							{
								auto t_id_out = get_task_id(m_ref_out.tasks, t_ref_out);
								auto s_id_out = get_socket_id(t_ref_out.sockets, *s_ref_out);

								assert(t_id_out != -1);
								assert(s_id_out != -1);

								(*this->all_modules[thread_id][m_id_out]).tasks[t_id_out]->set_autoalloc(true);

								auto &s_in  = *this->all_modules[thread_id][m_id    ]->tasks[t_id    ]->sockets[s_id    ];
								auto &s_out = *this->all_modules[thread_id][m_id_out]->tasks[t_id_out]->sockets[s_id_out];
								s_in.bind(s_out);
							}
						}
					}
				}
			}

			// add the sub-sequence to the current tree node
			sequence_cpy->set_contents(ss_cpy);

			for (size_t c = 0; c < sequence_ref->get_children().size(); c++)
				sequence_cpy->add_child(new Generic_node<Sub_sequence>(sequence_cpy,
				                                                       {},
				                                                       nullptr,
				                                                       sequence_cpy->get_depth() +1,
				                                                       0,
				                                                       c));

			for (size_t c = 0; c < sequence_ref->get_children().size(); c++)
				duplicate_sequence(sequence_ref->get_children()[c], sequence_cpy->get_children()[c], thread_id);
		}
	};

	std::function<void(Generic_node<Sub_sequence>*)> set_autoalloc_true =
		[&set_autoalloc_true](Generic_node<Sub_sequence>* node)
		{
			for (auto t : node->get_c()->tasks)
				t->set_autoalloc(true);
			for (auto c : node->get_children())
				set_autoalloc_true(c);
		};

	for (size_t thread_id = (this->tasks_inplace ? 1 : 0); thread_id < this->sequences.size(); thread_id++)
	{
		if (this->is_thread_pinning())
			Thread_pinning::pin(this->puids[thread_id]);

		this->sequences[thread_id] = new Generic_node<Sub_sequence>(nullptr, {}, nullptr, 0, 0, 0);
		duplicate_sequence(sequence, this->sequences[thread_id], thread_id);
		set_autoalloc_true(this->sequences[thread_id]);

		if (this->is_thread_pinning())
			Thread_pinning::unpin();
	}
}

template void tools::Sequence::duplicate<tools::Sub_sequence_const, const module::Module>(const Generic_node<tools::Sub_sequence_const>*);
template void tools::Sequence::duplicate<tools::Sub_sequence,             module::Module>(const Generic_node<tools::Sub_sequence      >*);

template <class SS>
void Sequence
::delete_tree(Generic_node<SS> *node)
{
	if (node != nullptr)
	{
		for (auto c : node->get_children())
			this->delete_tree(c);
		auto c = node->get_c();
		if (c != nullptr) delete c;
		delete node;
	}
}

template void tools::Sequence::delete_tree<tools::Sub_sequence_const>(Generic_node<tools::Sub_sequence_const>*);
template void tools::Sequence::delete_tree<tools::Sub_sequence      >(Generic_node<tools::Sub_sequence      >*);

template <class VTA>
void Sequence
::export_dot_subsequence(const VTA &subseq,
                         const std::vector<size_t> &tasks_id,
                         const subseq_t &subseq_type,
                         const std::string &subseq_name,
                         const std::string &tab,
                               std::ostream &stream) const
{
	if (!subseq_name.empty())
	{
		stream << tab << "subgraph \"cluster_" << subseq_name << "_" << +this << "\" {" << std::endl;
		stream << tab << tab << "node [style=filled];" << std::endl;
	}
	size_t exec_order = 0;
	for (auto &t : subseq)
	{
		std::string color = dynamic_cast<module::Adaptor*>(&t->get_module()) ? "green" :"blue";
		color = dynamic_cast<module::AProbe*>(&t->get_module()) ? "pink" : color;
		stream << tab << tab << "subgraph \"cluster_" << +&t->get_module() << "_" << +t << "\" {" << std::endl;
		stream << tab << tab << tab << "node [style=filled];" << std::endl;
		stream << tab << tab << tab << "subgraph \"cluster_" << +&t << "\" {" << std::endl;
		stream << tab << tab << tab << tab << "node [style=filled];" << std::endl;
		for (auto &s : t->sockets)
		{
			std::string stype;
			switch (t->get_socket_type(*s))
			{
				case module::socket_t::SIN: stype = "in"; break;
				case module::socket_t::SOUT: stype = "out"; break;
				case module::socket_t::SIN_SOUT: stype = "in_out"; break;
				default: stype = "unkn"; break;
			}
			stream << tab << tab << tab << tab << "\"" << +s.get() << "\""
			                                   << "[label=\"" << stype << ":" << s->get_name() << "\"];" << std::endl;
		}
		stream << tab << tab << tab << tab << "label=\"" << t->get_name() << " (id = " << tasks_id[exec_order] << ")" << "\";" << std::endl;
		stream << tab << tab << tab << tab << "color=" << color << ";" << std::endl;
		stream << tab << tab << tab << "}" << std::endl;
		stream << tab << tab << tab << "label=\"" << t->get_module().get_name() << "\n"
		                            << (t->get_module().get_custom_name().empty() ? "" : t->get_module().get_custom_name() + "\n")
		                            << "exec order: [" << exec_order++ << "]\n"
		                            << "addr: " << +&t->get_module() << "\";" << std::endl;
		stream << tab << tab << tab << "color=" << color << ";" << std::endl;
		stream << tab << tab << "}" << std::endl;
	}
	if (!subseq_name.empty())
	{
		stream << tab << tab << "label=\"" << subseq_name << "\";"
		       << std::endl;
		std::string color = subseq_type == subseq_t::LOOP ? "red" : "blue";
		stream << tab << tab << "color=" << color << ";" << std::endl;
		stream << tab << "}" << std::endl;
	}
}

template void tools::Sequence::export_dot_subsequence<std::vector<      module::Task*>>(const std::vector<      module::Task*>&, const std::vector<size_t>&, const subseq_t&, const std::string&, const std::string&, std::ostream&) const;
template void tools::Sequence::export_dot_subsequence<std::vector<const module::Task*>>(const std::vector<const module::Task*>&, const std::vector<size_t>&, const subseq_t&, const std::string&, const std::string&, std::ostream&) const;

template <class VTA>
void Sequence
::export_dot_connections(const VTA &subseq,
                         const std::string &tab,
                               std::ostream &stream) const
{
	for (auto &t : subseq)
	{
		for (auto &s : t->sockets)
		{
			if (t->get_socket_type(*s) == module::socket_t::SOUT ||
				t->get_socket_type(*s) == module::socket_t::SIN_SOUT)
			{
				for (auto &bs : s->get_bound_sockets())
				{
					stream << tab << "\"" << +s.get() << "\" -> \"" << +bs << "\"" << std::endl;
				}
			}
		}
	}
}

template void tools::Sequence::export_dot_connections<std::vector<      module::Task*>>(const std::vector<      module::Task*>&, const std::string&, std::ostream&) const;
template void tools::Sequence::export_dot_connections<std::vector<const module::Task*>>(const std::vector<const module::Task*>&, const std::string&, std::ostream&) const;

void Sequence
::export_dot(std::ostream &stream) const
{
	auto root = this->sequences[0];
	this->export_dot(root, stream);
}

template <class SS>
void Sequence
::export_dot(Generic_node<SS>* root, std::ostream &stream) const
{
	std::function<void(Generic_node<SS>*,
	                   const std::string&,
	                   std::ostream&)> export_dot_subsequences_recursive =
		[&export_dot_subsequences_recursive, this](Generic_node<SS>* cur_node,
		                                           const std::string &tab,
		                                           std::ostream &stream)
		{
			if (cur_node != nullptr)
			{
				this->export_dot_subsequence(cur_node->get_c()->tasks,
				                             cur_node->get_c()->tasks_id,
				                             cur_node->get_c()->type,
				                             "Sub-sequence "+std::to_string(cur_node->get_c()->id),
				                             tab,
				                             stream);

				for (auto c : cur_node->get_children())
					export_dot_subsequences_recursive(c, tab, stream);
			}
		};

	std::function<void(Generic_node<SS>*,
	                   const std::string&,
	                   std::ostream&)> export_dot_connections_recursive =
		[&export_dot_connections_recursive, this](Generic_node<SS> *cur_node,
		                                          const std::string &tab,
		                                          std::ostream &stream)
		{
			if (cur_node != nullptr)
			{
				this->export_dot_connections(cur_node->get_c()->tasks, tab, stream);

				for (auto c : cur_node->get_children())
					export_dot_connections_recursive(c, tab, stream);
			}
		};

	std::string tab = "\t";
	stream << "digraph Sequence {" << std::endl;
	export_dot_subsequences_recursive(root, tab, stream);
	export_dot_connections_recursive (root, tab, stream);
	stream << "}" << std::endl;
}

void Sequence
::gen_processes(const bool no_copy_mode)
{
	std::function<void(Generic_node<Sub_sequence>*)> gen_processes_recursive =
		[&gen_processes_recursive, no_copy_mode](Generic_node<Sub_sequence>* cur_node)
		{
			if (cur_node != nullptr)
			{
				std::map<module::Task*, std::function<int()>> modified_tasks;
				auto contents = cur_node->get_c();
				contents->processes.clear();
				contents->rebind_sockets.clear();
				contents->rebind_dataptrs.clear();
				for (auto task : contents->tasks)
				{
					if (dynamic_cast<module::Adaptor*>(&task->get_module()) &&
						task->get_name().find("pull") != std::string::npos && no_copy_mode)
					{
						auto pull_task = task;
						auto adp_pull = dynamic_cast<module::Adaptor*>(&pull_task->get_module());
						adp_pull->set_no_copy_pull(true);
						const auto rebind_id = contents->rebind_sockets.size();
						contents->rebind_sockets.resize(rebind_id +1);
						contents->rebind_dataptrs.resize(rebind_id +1);

						for (size_t s = 0; s < pull_task->sockets.size() -1; s++)
						{
							if (pull_task->get_socket_type(*pull_task->sockets[s]) == module::socket_t::SOUT)
							{
								std::vector<module::Socket*> bound_sockets;
								std::vector<void*> dataptrs;

								bound_sockets.push_back(pull_task->sockets[s].get());
								dataptrs.push_back(pull_task->sockets[s]->get_dataptr());

								auto bs = pull_task->sockets[s]->get_bound_sockets();
								bound_sockets.insert(bound_sockets.end(), bs.begin(), bs.end());
								for (auto sck : bs)
								{
									if (sck->get_task().get_socket_type(*sck) == module::socket_t::SIN_SOUT)
									{
										std::stringstream message;
										message << "'SIN_SOUT' socket type is unsupported at this time ("
										        << "'sck->get_name()' = " << sck->get_name() << ", "
										        << "'sck->get_task().get_name()' = ." << sck->get_task().get_name()
										        << ".)";
										throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
									}
									dataptrs.push_back(sck->get_dataptr());
								}

								contents->rebind_sockets[rebind_id].push_back(bound_sockets);
								contents->rebind_dataptrs[rebind_id].push_back(dataptrs);
							}
						}

						modified_tasks[pull_task] = [contents, pull_task, adp_pull, rebind_id]() -> int
						{
							// active or passive waiting here
							int ret = pull_task->exec();
							// rebind input sockets on the fly
							for (size_t sin_id = 0; sin_id < contents->rebind_sockets[rebind_id].size(); sin_id++)
							{
								if (contents->rebind_sockets[rebind_id][sin_id].size() > 1)
								{
									// we start to 1 because the rebinding of the 'pull_task' is made in the
									// 'pull_task->exec()' call (this way the debug mode is still working)
									auto swap_buff = contents->rebind_sockets[rebind_id][sin_id][1]->get_dataptr();
									auto buff = adp_pull->get_filled_buffer(sin_id, swap_buff);
									contents->rebind_sockets[rebind_id][sin_id][1]->bind(buff);
									// for the next tasks the same buffer 'buff' is required, an easy mistake is to re-swap
									// and the data will be false, this is why we just bind 'buff'
									for (size_t ta = 2; ta < contents->rebind_sockets[rebind_id][sin_id].size(); ta++)
										contents->rebind_sockets[rebind_id][sin_id][ta]->bind(buff);
								}
							}
							adp_pull->wake_up_pusher();
							return ret;
						};
					}

					if (dynamic_cast<module::Adaptor*>(&task->get_module()) &&
						task->get_name().find("push") != std::string::npos && no_copy_mode)
					{
						auto push_task = task;
						auto adp_push = dynamic_cast<module::Adaptor*>(&push_task->get_module());
						adp_push->set_no_copy_push(true);
						const auto rebind_id = contents->rebind_sockets.size();
						contents->rebind_sockets.resize(rebind_id +1);
						contents->rebind_dataptrs.resize(rebind_id +1);

						for (size_t s = 0; s < push_task->sockets.size() -1; s++)
							if (push_task->get_socket_type(*push_task->sockets[s]) == module::socket_t::SIN)
							{
								std::vector<module::Socket*> bound_sockets;
								std::vector<void*> dataptrs;

								bound_sockets.push_back(push_task->sockets[s].get());
								dataptrs.push_back(push_task->sockets[s]->get_dataptr());

								auto bound_socket = &push_task->sockets[s]->get_bound_socket();
								if (bound_socket->get_task().get_socket_type(*bound_socket) ==
								    module::socket_t::SIN_SOUT)
								{
									std::stringstream message;
									message << "'SIN_SOUT' socket type is unsupported at this time ("
									        << "'bound_socket->get_name()' = " << bound_socket->get_name() << ", "
									        << "'bound_socket->get_task().get_name()' = ."
									        << bound_socket->get_task().get_name() << ".)";
									throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
								}
								bound_sockets.push_back(bound_socket);
								dataptrs.push_back(bound_socket->get_dataptr());

								auto &bs = bound_socket->get_bound_sockets();
								for (size_t s = 0; s < bs.size(); s++)
									if (&bs[s]->get_task() != push_task)
									{
										if (bs[s]->get_task().get_socket_type(*bs[s]) == module::socket_t::SIN_SOUT)
										{
											std::stringstream message;
											message << "'SIN_SOUT' socket type is unsupported at this time ("
											        << "'s' = " << s << ", "
											        << "'bs[s]->get_name()' = " << bs[s]->get_name() << ", "
											        << "'bs[s]->get_task().get_name()' = ."
											        << bs[s]->get_task().get_name() << ".)";
											throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
										}

										bound_sockets.push_back(bs[s]);
										dataptrs.push_back(bs[s]->get_dataptr());
									}

								contents->rebind_sockets[rebind_id].push_back(bound_sockets);
								contents->rebind_dataptrs[rebind_id].push_back(dataptrs);
							}

						modified_tasks[push_task] = [contents, push_task, adp_push, rebind_id]() -> int
						{
							// active or passive waiting here
							auto ret = push_task->exec();
							// rebind output sockets on the fly
							for (size_t sout_id = 0; sout_id < contents->rebind_sockets[rebind_id].size(); sout_id++)
							{
								// we start to 1 because the rebinding of the 'push_task' is made in the
								// 'push_task->exec()' call (this way the debug mode is still working)
								auto swap_buff = contents->rebind_sockets[rebind_id][sout_id][1]->get_dataptr();
								auto buff = adp_push->get_empty_buffer(sout_id, swap_buff);
								contents->rebind_sockets[rebind_id][sout_id][1]->bind(buff);
								// the output socket linked to the push adp can have more than one socket bound and so
								// we have to rebind all the input sokects bound to the current output socket
								for (size_t ta = 2; ta < contents->rebind_sockets[rebind_id][sout_id].size(); ta++)
									contents->rebind_sockets[rebind_id][sout_id][ta]->bind(buff);
							}
							adp_push->wake_up_puller();
							return ret;
						};
					}
				}

				for (auto task : contents->tasks)
					if (modified_tasks.count(task))
						contents->processes.push_back(modified_tasks[task]);
					else
						contents->processes.push_back([task]() -> int { return task->exec(); });

				for (auto c : cur_node->get_children())
					gen_processes_recursive(c);
			}
		};

	size_t thread_id = 0;
	for (auto &sequence : this->sequences)
	{
		if (this->is_thread_pinning())
			Thread_pinning::pin(this->puids[thread_id++]);

		gen_processes_recursive(sequence);

		if (this->is_thread_pinning())
			Thread_pinning::unpin();
	}
}

void Sequence
::reset_no_copy_mode()
{
	std::function<void(Generic_node<Sub_sequence>*)> reset_no_copy_mode_recursive =
		[&reset_no_copy_mode_recursive](Generic_node<Sub_sequence>* cur_node)
		{
			if (cur_node != nullptr)
			{
				auto contents = cur_node->get_c();
				for (auto task : contents->tasks)
				{
					if (dynamic_cast<module::Adaptor*>(&task->get_module()) &&
						task->get_name().find("pull") != std::string::npos)
					{
						auto pull_task = task;
						auto adp_pull = dynamic_cast<module::Adaptor*>(&pull_task->get_module());
						adp_pull->set_no_copy_pull(false);
						adp_pull->reset_buffer();
					}

					if (dynamic_cast<module::Adaptor*>(&task->get_module()) &&
						task->get_name().find("push") != std::string::npos)
					{
						auto push_task = task;
						auto adp_push = dynamic_cast<module::Adaptor*>(&push_task->get_module());
						adp_push->set_no_copy_push(false);
						adp_push->reset_buffer();
					}
				}

				for (size_t rebind_id = 0; rebind_id < contents->rebind_sockets.size(); rebind_id++)
					for (size_t s = 0; s < contents->rebind_sockets[rebind_id].size(); s++)
						for (size_t ta = 0; ta < contents->rebind_sockets[rebind_id][s].size(); ta++)
							contents->rebind_sockets[rebind_id][s][ta]
								->bind(contents->rebind_dataptrs[rebind_id][s][ta]);
			}

			for (auto c : cur_node->get_children())
				reset_no_copy_mode_recursive(c);
		};

	for (auto &sequence : this->sequences)
		reset_no_copy_mode_recursive(sequence);

}

void Sequence
::set_no_copy_mode(const bool no_copy_mode)
{
	this->no_copy_mode = no_copy_mode;
}

bool Sequence
::is_no_copy_mode() const
{
	return this->no_copy_mode;
}

Sub_sequence* Sequence
::get_last_subsequence(const size_t tid)
{
	std::function<Sub_sequence*(Generic_node<Sub_sequence>*)> get_last_subsequence_recursive =
		[&get_last_subsequence_recursive](Generic_node<Sub_sequence>* cur_node) -> Sub_sequence*
		{
			if (cur_node != nullptr)
			{
				Sub_sequence* last_ss = nullptr;
				for (auto c : cur_node->get_children())
					last_ss = get_last_subsequence_recursive(c);
				return last_ss ? last_ss : cur_node->get_contents();
			}
			else
			{
				// should never happen
				return nullptr;
			}
		};

	return get_last_subsequence_recursive(this->sequences[tid]);
}

void Sequence
::update_tasks_id(const size_t tid)
{
	std::function<void(Generic_node<Sub_sequence>*, size_t&)> update_tasks_id_recursive =
		[&update_tasks_id_recursive](Generic_node<Sub_sequence>* cur_node, size_t& taid)
		{
			if (cur_node != nullptr)
			{
				Sub_sequence* ss = cur_node->get_contents();
				ss->tasks_id.resize(ss->tasks.size());
				std::iota(ss->tasks_id.begin(), ss->tasks_id.end(), taid);
				taid += ss->tasks_id.size();

				for (auto c : cur_node->get_children())
					update_tasks_id_recursive(c, taid);
			}
		};

	size_t taid = 0;
	return update_tasks_id_recursive(this->sequences[tid], taid);
}

std::vector<module::Task*> Sequence
::get_tasks_from_id(const size_t taid)
{
	std::function<void(Generic_node<Sub_sequence>*, const size_t, std::vector<module::Task*>&)> get_tasks_from_id_recursive =
		[&get_tasks_from_id_recursive](Generic_node<Sub_sequence>* cur_node,
		                               const size_t taid,
		                               std::vector<module::Task*> &tasks)
		{
			if (cur_node != nullptr)
			{
				Sub_sequence* ss = cur_node->get_contents();
				bool found = false;
				for (size_t t = 0; t < ss->tasks_id.size(); t++)
					if (ss->tasks_id[t] == taid)
					{
						tasks.push_back(ss->tasks[t]);
						found = true;
						break;
					}

				if (!found)
					for (auto c : cur_node->get_children())
						get_tasks_from_id_recursive(c, taid, tasks);
			}
		};

	std::vector<module::Task*> tasks;
	for (auto &s : this->sequences)
		get_tasks_from_id_recursive(s, taid, tasks);
	return tasks;
}

void Sequence
::update_firsts_and_lasts_tasks()
{
	this->firsts_tasks.clear();
	this->firsts_tasks.resize(this->n_threads);
	for (auto taid : firsts_tasks_id)
	{
		auto tasks = this->get_tasks_from_id(taid);
		for (size_t tid = 0; tid < tasks.size(); tid++)
			firsts_tasks[tid].push_back(tasks[tid]);
	}

	this->lasts_tasks.clear();
	this->lasts_tasks.resize(this->n_threads);
	for (auto taid : lasts_tasks_id)
	{
		auto tasks = this->get_tasks_from_id(taid);
		for (size_t tid = 0; tid < tasks.size(); tid++)
			lasts_tasks[tid].push_back(tasks[tid]);
	}
}