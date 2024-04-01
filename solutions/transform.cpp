#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

int main() {

    std::vector<double> v{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

    std::transform(v.begin(), v.end(), v.begin(), [] (double e) {return std::sqrt(e);});

    for (auto e : v) {
	std::cout << e << std::endl;
    }
    
    auto sum = std::accumulate(v.cbegin(), v.cend(), 0.0);

    std::cout << sum << std::endl;
}
