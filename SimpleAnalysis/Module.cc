#include "Module.h"

Module::Module(ExRootTreeReader* data)
{
  _data = data;
}

Module::~Module()
{
}

void Module::initialize() 
{
}

bool Module::execute(std::map<std::string, std::any>* DataStore)
{
  return true;
}

void Module::finalize()
{
}


