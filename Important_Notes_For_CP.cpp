/********   All Required Header Files ********/

#include <bits/stdc++.h> // This Contain All Header Files

using namespace std;

/****** Fast I/O Methods *********/

void fastio() {
	ios_base::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);
}

/************************************/


/******** User-defined Function *******/

void solve() {

}

/**************************************/


/********** Main()  function **********/
int main()
{
	//freopen("input.txt", "r", stdin);
	//freopen("output.txt","w",stdout);

	int tc;
	cin >> tc;

	while (tc--) 
		solve();

	return 0;
}
/********  Main() Ends Here *************/

/************ Important Function **************/ 

/*

	1. `equal()`
		
		Use: Checks if two ranges are equal.

    2. `iota()`
		
		Use: Fills a range with sequentially increasing values.

	3. `max_element()` & `min_element()`
		
		Use: Finds the largest/smallest element in a range.

	4. `distance()`
		
		Use: Calculates the number of elements between two iterators.

	5. `count()`
		
		Use: Counts the number of elements that are equal to a given value.

	6. `replace()` & `replace_if()`
		
		Use: Replaces values in a range.

	7. `search()`
		
		Use: Searches for a subsequence in a range.

	8. `for_each()`
		
		Use: Applies a function to each element in a range.

	9. `sort()`
		
		Use: Sorts elements in a range.

	10. `accumulate()`
		
		Use: Calculates the sum of elements in a range.

	11. `inner_product()`
		
		Use: Calculates the inner product of two ranges.

	12. `partial_sum()`
		
		Use: Calculates the partial sums of elements in a range.

	13. `adjacent_difference()` & `adjacent_find()`
		
		Use: Computes the differences between adjacent elements / finds first occurrence of adjacent equal elements.

	 14. `rotate()`
		
		Use: Rotates elements in a range to the left.

	15. `shuffle()`
		
		Use: Randomly shuffles elements in a range.

	16. `next_permutation()`
		
		Use: Generates the next lexicographical permutation.

	17. `binary_search()`
		
		Use: Checks if a value exists in a sorted range.

	18. `lower_bound()` & `upper_bound()`
		
		Use: Finds the first/last position where a value could be inserted to maintain order.

	19. `unique()`
		
		Use: Removes consecutive duplicates in a range.

	20. `set_intersection()` & `set_union()`
		
		Use: Computes the intersection/union of two sorted ranges.

	21. `set_difference()`
		
		Use: Computes the difference between two sorted ranges.

	22. `set_symmetric_difference()`
		
		Use: Computes the symmetric difference between two sorted ranges.

*/

int main() {

	// Example usage of equal()
	std::vector<int> a = { 1, 2, 3 };
	std::vector<int> b = { 1, 2, 3 };
	bool result = std::equal(a.begin(), a.end(), b.begin()); // result is true

	// Example usage of iota()
	std::vector<int> v(5);
	std::iota(v.begin(), v.end(), 10); // v is {10, 11, 12, 13, 14}

	// Example usage of max_element() & min_element()
	std::vector<int> v = { 3, 1, 4, 1, 5 };
	auto max_it = std::max_element(v.begin(), v.end()); // *max_it is 5
	auto min_it = std::min_element(v.begin(), v.end()); // *min_it is 1

	// Example usage of distance()
	std::vector<int> v = { 1, 2, 3, 4, 5 };
	auto dist = std::distance(v.begin(), v.end()); // dist is 5

	// Example usage of count()
	std::vector<int> v = { 1, 2, 3, 2, 4 };
	int count = std::count(v.begin(), v.end(), 2); // count is 2

	// Example usage of replace() & replace_if()
	std::vector<int> v = { 1, 2, 3, 2, 4 };
	std::replace(v.begin(), v.end(), 2, 5); // v is {1, 5, 3, 5, 4}
	std::replace_if(v.begin(), v.end(), [](int n) { return n > 3; }, 0); // v is {1, 0, 3, 0, 0}

	// Example usage of search()
	std::vector<int> v = { 1, 2, 3, 4, 5 };
	std::vector<int> sub = { 3, 4 };
	auto it = std::search(v.begin(), v.end(), sub.begin(), sub.end()); // it points to 3 in v

	// Example usage of for_each()
	std::vector<int> v = { 1, 2, 3, 4, 5 };
	std::for_each(v.begin(), v.end(), [](int& n) { n *= 2; }); // v is {2, 4, 6, 8, 10}

	// Example usage of sort()
	std::vector<int> v = { 3, 1, 4, 1, 5 };
	std::sort(v.begin(), v.end()); // v is {1, 1, 3, 4, 5}

	// Example usage of accumulate()
	std::vector<int> v = { 1, 2, 3, 4, 5 };
	int sum = std::accumulate(v.begin(), v.end(), 0); // sum is 15

	// Example usage of inner_product()
	std::vector<int> a = { 1, 2, 3 };
	std::vector<int> b = { 4, 5, 6 };
	int product = std::inner_product(a.begin(), a.end(), b.begin(), 0); // product is 32

	// Example usage of partial_sum()
	std::vector<int> v = { 1, 2, 3, 4 };
	std::vector<int> result(4);
	std::partial_sum(v.begin(), v.end(), 0); // result is {1, 3, 6, 10}

	// Example usage of adjacent_difference() & adjacent_find()
	std::vector<int> v = { 1, 3, 6, 10 };
	std::vector<int> diff(4);
	std::adjacent_difference(v.begin(), v.end(), diff.begin()); // diff is {1, 2, 3, 4}
	auto it = std::adjacent_find(v.begin(), v.end()); // it points to end (no adjacent elements)

	// Example usage of rotate()
	std::vector<int> v = { 1, 2, 3, 4, 5 };
	std::rotate(v.begin(), v.begin() + 2, v.end()); // v is {3, 4, 5, 1, 2}

	// Example usage of shuffle()
	std::vector<int> v = { 1, 2, 3, 4, 5 };
	//std::random_device rd;
	void fun();
	std::shuffle(v.begin(), v.end(), fun); // v is shuffled randomly

	// Example usage of next_permutation()
	std::vector<int> v = { 1, 2, 3 };
	std::next_permutation(v.begin(), v.end()); // v is {1, 3, 2}

	// Example usage of binary_search()
	std::vector<int> v = { 1, 2, 3, 4, 5 };
	bool found = std::binary_search(v.begin(), v.end(), 3); // found is true

	// Example usage of lower_bound() & upper_bound()
	std::vector<int> v = { 1, 2, 2, 3, 4 };
	auto lb = std::lower_bound(v.begin(), v.end(), 2); // lb points to the first 2
	auto ub = std::upper_bound(v.begin(), v.end(), 2); // ub points to the element after the last 2

	// Example usage of unique()
	std::vector<int> v = { 1, 2, 2, 3, 3, 4 };
	auto it = std::unique(v.begin(), v.end()); // v is {1, 2, 3, 4, ?, ?}
	v.erase(it, v.end()); // v is {1, 2, 3, 4}

	// Example usage of set_intersection() & set_union()
	std::vector<int> a = { 1, 2, 3 };
	std::vector<int> b = { 2, 3, 4 };
	std::vector<int> intersect(2);
	std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), intersect.begin()); // intersect is {2, 3}
	std::vector<int> uni(4);
	std::set_union(a.begin(), a.end(), b.begin(), b.end(), uni.begin()); // uni is {1, 2, 3, 4}

	// Example usage of set_difference()
	std::vector<int> a = { 1, 2, 3 };
	std::vector<int> b = { 2, 3, 4 };
	std::vector<int> diff(1);
	std::set_difference(a.begin(), a.end(), b.begin(), b.end(), diff.begin()); // diff is {1}

	// Example usage of set_symmetric_difference()
	std::vector<int> a = { 1, 2, 3 };
	std::vector<int> b = { 2, 3, 4 };
	std::vector<int> sym_diff(2);
	std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(), sym_diff.begin()); // sym_diff is {1, 4}
}
