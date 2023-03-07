#pragma once
#ifndef LI_Decay_H
#define LI_Decay_H

namespace LI {
namespace processes {

class Decay : public Process{
  friend cereal::access;
  private:
  public: 

}; // class Decay

} // namespace crosssections
} // namespace LI

CEREAL_CLASS_VERSION(LI::processes::Decay, 0);

#endif // LI_Decay_H
