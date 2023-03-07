#pragma once
#ifndef LI_Decay_H
#define LI_Decay_H

namespace LI {
namespace processes {

class Decay{
  friend cereal::access;
  private:
  public: 
    Decay();
    virtual ~Decay() {};

}; // class Decay

} // namespace processes
} // namespace LI

CEREAL_CLASS_VERSION(LI::processes::Decay, 0);

#endif // LI_Decay_H
