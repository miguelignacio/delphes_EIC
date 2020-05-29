#ifndef ANALYSISFUNCTIONS
#define ANALYSISFUNCTIONS

#include <vector>

template <class T> std::vector<T*> SelectorFcn(std::vector<T*> particles, bool selector(T*)) 
{
  std::vector<T*> output_list;
  for (T* particle : particles) {
    if (selector(particle)) {
      output_list.push_back(particle);
    }
  }

  return output_list;
}


#endif
