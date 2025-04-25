
%module bishop
%{
#include "NoiseMachine.h"
%}

//%include <std_shared_ptr.i>
%include <boost_shared_ptr.i>
%shared_ptr(NoiseMachineBase);

%include "NoiseMachine.h"

%template(NoiseMachineBase)  std::shared_ptr<lux::Noise>;
