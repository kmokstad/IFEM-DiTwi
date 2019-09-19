#include <algorithm>
#include <array>
#include <map>
#include <numeric>
#include <vector>

template<class Key, class Value>
struct NelderMead
{
  struct Entry {
    Key key;
    Value val;
  };

  template<class Function, class Convergence>
  static Key optimize(std::vector<Entry>& vectors, Value& result, int dim, double tol,
                      Function&& f, Convergence&& conv)
  {
    bool done = conv(vectors[0].val) < tol;

    while (!done) {
      std::sort(vectors.begin(), vectors.end(), [f](const Entry& a, const Entry& b)
                                                {
                                                  return a.val < b.val;
                                                });
      Key cog(vectors[0].key);
      cog.fill(0.0);
      for (size_t i = 0; i < vectors.size()-1; ++i)
        cog += vectors[i].key;

      cog /= vectors.size()-1;

      auto best = vectors.front();
      auto worst = vectors.back();
      auto second_worst = vectors[vectors.size()-2];

      // reflect
      const double alpha = 1.0;
      const double gamma = 2.0;
      const double rho   = 0.5;
      const double sigma = 0.5;

      auto reflected = cog + alpha*(cog - worst.key);
      auto fr = f(reflected);
      const auto& fb = best.val;
      const auto& fsw = second_worst.val;

      if (fr < fsw && fr > fb) {
        vectors.back().key = reflected;
        vectors.back().val = fr;
      } else if (fr < fb) {
        // expand
        auto expanded = cog + gamma*(reflected - cog);
        auto fe = f(expanded);
        if (fe < fr) {
          vectors.back().key = expanded;
          vectors.back().val = fe;
        } else {
          vectors.back().key = reflected;
          vectors.back().val = fr;
        }
      } else {
        // contract
        auto contracted = cog + rho*(worst.key - cog);
        auto fc = f(contracted);
        if (fc < worst.val) {
          vectors.back().key = contracted;
          vectors.back().val = fc;
        } else {
          for (size_t i = 0; i < vectors.size(); ++i) {
            vectors[i].key = best.key + sigma*(vectors[i].key - best.key);
            vectors[i].val = f(vectors[i].key);
          }
        }
      }
      done = conv(fb) < tol;
      if (done)
        result = fb;
    }

    return vectors[0].key;
  } 
};
