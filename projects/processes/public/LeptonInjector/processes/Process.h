#pragma once
#ifndef LI_Process_H
#define LI_Process_H

namespace LI {
namespace processes {

class Process {
  friend cereal::access;
  private:
  public: 

}; // class Process

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::processes::Process, 0);

#endif // LI_Process_H
