/*    
	Copyright (C) 2021
	This minimal hash table library was written in C++ based on the standard Python(2.7.6) library for a dictionary object:
		https://www.python.org/downloads/release/python-276/
		All Python derivative work was written based on Gnu Public License compatible material licensed to the Python Software Foundation.
	All rights reserved.
	
	License: BSD (Berkeley Software Distribution)
	
	Author: Chris Burgoyne <christopher.burgoyne@umontana.edu>
	***Base class for dicitonary object used in DiGraph.h, a DiGraph library written for the SweetSEQer algorithm in C++***
*/

#ifndef Dict_H
#define Dict_H

#include<cmath>
#include<cassert>

using namespace std;

template<typename item_key, typename item_value>
struct dict_item{
	int used = 0;
	item_key first = 0.0;
	item_value second;
};

#define MAX_VAL 2147483648.0
#define MINSIZE 8
#define PERTURB_SHIFT 5
#define PERTURB_A 2
#define PERTURB_B 1
#define RESIZEFACTOR 4
#define MAXFILL 2.0/3.0
#define DBL_MAX 1E+37
#define USED 1
#define EMPTY 0

template<typename dictKey, typename dictEntry>
class Dict {
	public:
	dict_item<dictKey, dictEntry> *table = new dict_item<dictKey, dictEntry>[MINSIZE];
	unsigned long table_size = MINSIZE, table_used = 0;
	int table_first_entry = MINSIZE, table_last_entry = 0;
	
	unsigned long hash_double(double v) {
		double intpart, fractpart;
		int expo;
		unsigned long hipart;
		unsigned long hash_value;
		fractpart = modf(v, &intpart);
		if (fractpart == 0.0) {
			hash_value = (long)intpart;
			return hash_value;
		}
		v = frexp(v, &expo);
		v *= MAX_VAL;
		hipart = (long)v;
		v = (v - (double)hipart) * MAX_VAL;
		hash_value = hipart + (long)v + (expo << 15);
		return hash_value;
	}
	
	unsigned long hash_double_pair(double u, double v) {
		long x, y;
		ssize_t len = 1;
		long mult = 1000003L;
		x = 0x345678L;
		
		y = hash_double(u);
		x = (x ^ y) * mult;
		mult += (long)(82520L + len + len);
		--len;
		y = hash_double(v);
		x = (x ^ y) * mult;
		mult += (long)(82520L + len + len);
		
		x += 97531L;
		if (x == -1)
			x = -2;
		return x;
	}
	
	void track_table_start_end(unsigned long insertion_idx) {
		// maintaining the iteration list for iteration in time proportional to (last used bin - first used bin)
		if (table_used == 1) {
			table_first_entry = insertion_idx;
			table_last_entry = insertion_idx;	
		}
		// new last entry
		else if ((int)insertion_idx > table_last_entry)
			table_last_entry = (int)insertion_idx;	
		// new first entry
		else if ((int)insertion_idx < table_first_entry)
			table_first_entry = (int)insertion_idx;
	}
	
	unsigned long find_insertion_idx(dictKey key){
		unsigned long hash = hash_double(key);
		unsigned long mask = table_size - 1;
		unsigned long i = hash & mask;
		unsigned long freeslot = -1;
		if (table[i & mask].first == key)
			return i & mask;
		if (table[i & mask].used == EMPTY)
			freeslot = i & mask;
		
		for (unsigned long perturb = hash; ; perturb >>= PERTURB_SHIFT) {
			i = (i << 2) + i + perturb + 1;
			if (table[i & mask].used == EMPTY && table[i & mask].first != 0.0)
				continue;
			if (table[i & mask].used == EMPTY)
				return freeslot == -1 ? (i & mask) : freeslot;
			if (table[i & mask].first == key)
				return i & mask;
		}
		// failed to find empty slot
		return -1;
	}
	
	unsigned long find_idx(dictKey key){
		unsigned long hash = hash_double(key);
		unsigned long mask = table_size - 1;
		unsigned long i = hash & mask;
		unsigned long perturb = hash;
		int counter = 0;
		
		if (table[i & mask].first == key)
			return i & mask;
		for ( ; counter < table_size; perturb >>= PERTURB_SHIFT) {
			i = (i * 5) + perturb + 1;
			if (table[i & mask].used == USED && table[i & mask].first == key) {
				return i & mask;
			}
			++counter;
		}
		return i & mask;
	}
	
	void insert_key_val_pair(dictKey key, dictEntry val) {
		// if key is new find an available location in the table, else find the existing location of the key
		unsigned long insertion_idx = find_insertion_idx(key);
		// new entry
		if (table[insertion_idx].used != USED) {
			// insert key-value pair into the table
			table[insertion_idx].used = USED;
			table[insertion_idx].first = key;
			table[insertion_idx].second = val;
			table_used++;
			// reset tracking of table first and last entries for iterator
			track_table_start_end(insertion_idx);
		// update existing entry
		} else
			table[insertion_idx].second = val;
	}
	
	void insert_new_dict_by_key(dictKey key) {
		unsigned long insertion_idx = find_insertion_idx(key);
		table[insertion_idx].used = USED;
		table[insertion_idx].first = key;
		table_used++;
		track_table_start_end(insertion_idx);
	}
	
	void insert_clean_dict(dictKey key, dictEntry *val) {
		// if key is new find an available location in the table, else find the existing location of the key
		unsigned long insertion_idx = find_insertion_idx(key);
		// insert key and deep copy existing dictionary contents
		table[insertion_idx].used = USED;
		table[insertion_idx].first = key;
		table[insertion_idx].second.copy_nodes(val);
		table_used++;
		// reset tracking of table first and last entries for iterator	
		track_table_start_end(insertion_idx);
	}
	
	void resize(ssize_t _newsize, dictKey multiplier) {
		int old_table_used = table_used, old_table_last = table_last_entry, old_table_first = table_first_entry;
		unsigned long old_size = table_size;
		table_used = 0;
		// find new table size using bit shifting where new_size is a power of 2 greater than (old_size * RESIZEFACTOR)
		for (table_size = MINSIZE; table_size < _newsize && table_size > 0; table_size <<= PERTURB_B)
			;
		// reset tracking for table first and last entries to relative initial values
		table_first_entry = table_size - 1, table_last_entry = 0;
		
		// save pointer to old data
		dict_item<dictKey, dictEntry> *old_table = table;
		
		// allocate new memory for table
		table = new dict_item<dictKey, dictEntry>[table_size];
		
		// insert all dict_entries into new memory
		for (int i = old_table_first; i <= old_table_last; ++i)
			if (old_table[i].used == USED)
				insert_key_val_pair(old_table[i].first, old_table[i].second);
		
		// clear the memory in the old table
		delete[] old_table;
	}
	
	void resize_nested_dictionary(ssize_t _newsize, dictKey multiplier) {
		int old_table_used = table_used, old_table_last = table_last_entry, old_table_first = table_first_entry;
		unsigned long old_size = table_size;
		table_used = 0;
		// find new table size using bit shifting where new_size is a power of 2 greater than (old_size * RESIZEFACTOR)
		for (table_size = MINSIZE; table_size < _newsize && table_size > 0; table_size <<= 1)
			;
		// reset tracking for table first and last entries to relative initial values
		table_first_entry = table_size - 1, table_last_entry = 0;
		
		// save pointer to old data
		dict_item<dictKey, dictEntry> *old_table = table;
		
		// allocate new memory for table
		table = new dict_item<dictKey, dictEntry>[table_size];
		
		// insert all dict_entries into new memory
		for (int i = old_table_first; i <= old_table_last; ++i)
			if (old_table[i].used == USED)
				insert_clean_dict(old_table[i].first, &old_table[i].second);

		// clear the memory in the old table
		delete[] old_table;
	}
	
	unsigned long find_insertion_idx_dbl_pair(dictKey key, double val){
		unsigned long hash = hash_double_pair(key,val);
		unsigned long mask = table_size - 1;
		unsigned long i = hash & mask;
		unsigned long perturb = hash;
		
		if (table[i & mask].used == EMPTY || (table[i & mask].first == key && table[i & mask].second == val))
			return i & mask;
		for ( ; ; perturb >>= PERTURB_SHIFT) {
			i = (i << 2) + i + perturb + 1;
			if (table[i & mask].used == EMPTY)
				return i & mask;
			if (table[i & mask].used == USED && table[i & mask].first == key && table[i & mask].second == val)
				return i & mask;
		}
		return i & mask;
	}
	
	void insert_pair_doubles(double key, double val) {
		unsigned long insertion_idx = find_insertion_idx_dbl_pair(key, val);
		
		if (table[insertion_idx].used != USED) {
			table[insertion_idx].used = USED;
			table[insertion_idx].first = key;
			table[insertion_idx].second = val;
			table_used++;
			
			track_table_start_end(insertion_idx);
		} else 
			table[insertion_idx].second = val;
	}
	
	void insert_clean_double_pair(double key, double val) {
		unsigned long insertion_idx = find_insertion_idx_dbl_pair(key, val);
		
		table[insertion_idx].used = USED;
		table[insertion_idx].first = key;
		table[insertion_idx].second = val;
		++table_used;
		
		track_table_start_end(insertion_idx);
	}
	
	void resize_double(ssize_t _newsize, dictKey multiplier) {
		ssize_t newsize;
		int old_table_used = table_used;
		table_used = 0;
		
		// calculate new size
		for (newsize = MINSIZE; newsize < _newsize && newsize > 0; newsize <<= 1)
			;
		unsigned long old_size = table_size;
		table_size = newsize;
		table_first_entry = newsize - 1, table_last_entry = 0;
		
		// save pointer to old data
		dict_item<dictKey, dictEntry> *old_table = table;
		
		// allocate new memory for table
		table = new dict_item<dictKey, dictEntry>[newsize];
		
		// insert all dict_entries into new memory
		for (int i = 0; i < old_size; ++i)
			if (old_table[i].used == USED)
				insert_clean_double_pair(old_table[i].first, old_table[i].second);

		// clear the memory in the old table
		delete[] old_table;
	}
	
	void _remove(dictKey key){
		// find position of key in the table
		int insertion_idx = -1;
		for (int i = table_first_entry; i <= table_last_entry; ++i)
			if (table[i].used == USED && table[i].first == key){
				insertion_idx = i;
				break;
			}
		// key not found
		if (insertion_idx == -1)
			return;
		// removal process
		table[insertion_idx].used = EMPTY;
		table_used--;
		// reset position of first or last for the iterator
		if (table_used == 0)
			table_first_entry = table_size, table_last_entry = 0;
		else {
			if ((int)insertion_idx == table_last_entry)
				while (table[table_last_entry].used != USED && table_last_entry > 0)
					--table_last_entry;
			if ((int)insertion_idx == table_first_entry)
				while (table[table_first_entry].used != USED && table_first_entry < table_size)
					++table_first_entry;
		}
	}
	
	public:
		Dict() {}
		
		~Dict() {
			delete[] table;	
		}

		void insert(dictKey key, dictEntry value) {
			insert_key_val_pair(key, value);
			if(table_used > table_size * MAXFILL) {
				resize(table_size * RESIZEFACTOR, RESIZEFACTOR);
			}
		}
		
		void insert(double key, double value, bool is_double) {
			insert_pair_doubles(key, value);
			if(table_used > table_size * MAXFILL)
				resize_double(table_size * RESIZEFACTOR, RESIZEFACTOR);
		}
		
		void insert(dictKey key) {
			insert_new_dict_by_key(key);
			if(table_used > table_size * MAXFILL) {
				resize_nested_dictionary(table_size * RESIZEFACTOR, RESIZEFACTOR);
			}
		}
		
		void update(Dict *from) {
				if (from->table_used == 0)
					return;
				if ((table_used + from->table_used)*3 >= table_size*2)
						resize((table_used + from->table_used) * 2, 2);
				for (auto& nodes : *from)
						insert(nodes.first, nodes.second);

				/* 
				ssize_t updateSize = from->size() + size(), newTableSize;
				if (updateSize > (int)(table_size * MAXFILL)) {
					ssize_t sizeOverUpdate, newSize;
					for (sizeOverUpdate = MINSIZE; sizeOverUpdate < updateSize; sizeOverUpdate <<= 1)
						;
					sizeOverUpdate <<= 1;
					for (newTableSize = MINSIZE; newTableSize < sizeOverUpdate; newTableSize <<= 1)
						;
					//cout << sizeOverUpdate << " " << newTableSize << "\n";
					if (table_used > 0) {
						dict_item<dictKey, dictEntry> *old_table = table;
						table = new dict_item<dictKey, dictEntry>[newTableSize];
						int old_table_first = table_first_entry, old_table_last = table_last_entry, old_table_size = table_size;
						table_used = 0, table_first_entry = newTableSize, table_last_entry = 0, table_size = newTableSize;
						if (((updateSize & (updateSize-1)) == 0) && updateSize != 0) {
							resize(table_size * 2, 2);
						}
						for (int i = old_table_first; i <= old_table_last; ++i)
							if (old_table[i].used == USED) {
								insert(old_table[i].first, old_table[i].second);
							}
						delete[] old_table;
					}
					for (auto& nodes : *from) {
						insert(nodes.first, nodes.second);
					}
				} else {
					if (((updateSize & (updateSize-1)) == 0) && updateSize != 0) {
						resize(table_size * 2, 2);
					}
					for (auto& nodes : *from)
						insert(nodes.first, nodes.second);
				} */
		}
		
		void copy_nodes(Dict *from) {
			// if table to copy is larger than MINSIZE(8), create new table of incoming table size
			if (from->table_size > table_size) {
				delete[] table;
				table_size = from->table_size;
				table = new dict_item<dictKey, dictEntry>[table_size];
			}
			table_used = from->table_used;
			// copy values
			for (int i = 0; i < table_size; ++i) {
				if (from->table[i].used == USED) {
					table[i].used = USED;
					table[i].first = from->table[i].first;
					table[i].second = from->table[i].second;
				}
			}
			
			// track first and last entries for iterator
			for (int i = 0; i < table_size; ++i)
				if (table[i].used == USED) {
					table_first_entry = i;
					break;
				}
			for (int i = table_size - 1; i > 0; --i)
				if (table[i].used == USED) {
					table_last_entry = i;
					break;
				}
		}
		
		void copy_nested_dict(Dict *from) {
			// if table to copy is larger than MINSIZE(8), create new table of incoming table size
			if (from->table_size > table_size) {
				delete[] table;
				table_size = from->table_size;
				table = new dict_item<dictKey, dictEntry>[table_size];
			}
			table_used = from->table_used;
			// copy values
			for (int i = 0; i < table_size; ++i) {
				if (from->table[i].used == USED) {
					table[i].used = USED;
					table[i].first = from->table[i].first;
					table[i].second.copy_nodes(&from->table[i].second);
				}
			}
			
			// track first and last entries for iterator
			for (int i = 0; i < table_size; ++i)
				if (table[i].used == USED) {
					table_first_entry = i;
					break;
				}
			for (int i = table_size - 1; i > 0; --i)
				if (table[i].used == USED) {
					table_last_entry = i;
					break;
				}
		}
		
		void erase(dictKey key) {
			_remove(key);
		}
		
		bool contains(dictKey key) {
			if (table_used > 0) {
				//unsigned long query_idx = find_idx(key);
				for (int query_idx = table_first_entry; query_idx <= table_last_entry; ++query_idx)
					if (table[query_idx].used == USED && table[query_idx].first == key)
						return &table[query_idx];
			}
			return false;
		}
		
		bool contains_pair(double key, double val) {
			if (table_used > 0) {
				unsigned long query_idx = find_insertion_idx_dbl_pair(key,val);
				if (table[query_idx].used == USED && table[query_idx].first == key && table[query_idx].second == val)
					return true;
			}
			return false;
		}
		
		bool contains_either_of_pair(double key, double val) {
			if (table_used > 0) {
				for (int i = table_first_entry; i <= table_last_entry; ++i)
					if (table[i].used == USED && (table[i].first == key || table[i].first == val))
						return true;
			}
			return false;
		}
		
		dict_item<dictKey, dictEntry>* get(dictKey key) {
			if (table_used > 0) {
				//unsigned long query_idx = find_idx(key);
				for (int query_idx = table_first_entry; query_idx <= table_last_entry; ++query_idx)
					if (table[query_idx].used == USED && table[query_idx].first == key)
						return &table[query_idx];
			}
			return NULL;
		}
		
		int size(){
			return table_used;
		}
		
		double min_entry() {
			double min = DBL_MAX;
			for (int i = table_first_entry; i <= table_last_entry; ++i)
				if (table[i].used == USED && table[i].first < min)
					min = table[i].first;
			return min;
		}
		
		double max_entry() {
			double max = 0;
			for (int i = table_first_entry; i <= table_last_entry; ++i)
				if (table[i].used == USED && table[i].first > max)
					max = table[i].first;
			return max;
		}
		
		void print(){
			for (int i = 0; i < table_size; ++i)
			if (table[i].used == USED)
				cout << table[i].first << " ";
			else cout << i << " ";
			cout << " end\n";
		}

		void print_used(){
			for (int i = 0; i < table_size; ++i)
				cout << i << " ";
			cout << "\n";
			for (int i = 0; i < table_size; ++i)
				cout << table[i].used << " ";
			cout << " end\n";
		}
		
		void print_all(){
			for (int i = 0; i < table_size; ++i)
				if (table[i].used == USED)
					cout << table[i].first << ":" << table[i].second << "  ";
				else cout << i << "  ";
			cout << " end\n";
		}
		
		void print_mem() {
			for (int i = 0; i < table_size; ++i)
				cout << "\t*" << i << "\t";
			cout << "\n";
			for (int i = 0; i < table_size; ++i)
					cout << &table[i] << "\t";
			cout << "\n";
		}
		
		struct Iterator {
			using pointer = dict_item<dictKey,dictEntry>*;
			using reference = dict_item<dictKey,dictEntry>&;
			public:
				dict_item<dictKey,dictEntry>* ptr;
				int table_sz, iterator;
				friend class Dict;
				dict_item<dictKey, dictEntry> *t;
				Iterator(pointer n, int m, int o, dict_item<dictKey, dictEntry> *tt=NULL) : ptr(n), iterator(m), table_sz(o), t(tt) {}
			
				Iterator operator++() {
					++ptr;
					++iterator;
					while (iterator < table_sz && ptr->used == false) {	
						++ptr;
						++iterator;
					}
					return *this;
				}
				friend bool operator!=(const Iterator& a, const Iterator& b) {
					return a.ptr != b.ptr;
				}
				reference operator*() { return *ptr; }
				pointer operator->() { return ptr; }
		};
		
		Iterator begin() {
			return Iterator(&table[table_first_entry], table_first_entry, table_size, table);
		}
		
		Iterator end()   {
			return Iterator(&table[table_size], table_first_entry, table_size);
		}
		
};


#endif
