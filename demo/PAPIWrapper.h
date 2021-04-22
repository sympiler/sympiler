//
// Created by kazem on 2020-05-26.
//

#ifndef FUSION_PAPIWRAPPER_H
#define FUSION_PAPIWRAPPER_H

#include <papi.h>
#include <utility>

namespace sym_lib{

 struct EventCounterBundle{
  int event_code;
  std::string event_name;
  std::vector<long long> counters;
 };

 class PAPIWrapper {
  int num_events_;
  int *event_list_; // event id list passed to PAPI
  long long *event_counters_; // output counter from PAPI
  int event_set_; // PAPI

  std::string convert_code_to_string(int code){
   char EventCodeStr[PAPI_MAX_STR_LEN];
   PAPI_event_code_to_name(code,EventCodeStr);
   std::string ret(EventCodeStr);
   return ret;
  }

 public:
  std::vector<EventCounterBundle> counter_bundles_;
  PAPIWrapper(const std::vector<int>& event_list, int num_instances){
   num_events_ = event_list.size();
   event_list_ = new int[num_events_]();
   event_counters_ = new long long [num_events_]();
   for (int i = 0; i < num_events_; ++i) {
    event_list_[i] = event_list[i];
    EventCounterBundle ecb;
    ecb.event_code = event_list_[i];
    ecb.event_name = convert_code_to_string(event_list_[i]);
    counter_bundles_.emplace_back(ecb);
   }
   event_set_ = PAPI_NULL;
  }

  ~PAPIWrapper(){
   delete []event_list_;
   delete []event_counters_;
  }

  int begin_profiling(){
   event_set_ = PAPI_NULL;
   int retval;
   if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
    std::cerr << "PAPI error initializing " << retval << "\n";
    return 0;
   }
   if((retval = PAPI_create_eventset(&event_set_)) != PAPI_OK) {
    std::cerr << "PAPI error creating event set " << retval << "\n";
    return 0;
   }
   if((retval = PAPI_add_events(event_set_, event_list_, num_events_)) != PAPI_OK) {
    std::cerr << "PAPI error adding events " << retval << "\n";
    return 0;
   }
   for (int i = 0; i < num_events_; ++i)
    event_counters_[i] = 0;
   if((retval = PAPI_reset(event_set_)) != PAPI_OK)
    std::cerr << "PAPI error " << retval << "\n";
   if((retval = PAPI_start(event_set_)) != PAPI_OK)
    std::cerr << "PAPI error " << retval << "\n";
  }

  void finish_profiling(){
   int retval=0;
   if((retval = PAPI_stop(event_set_, event_counters_)) != PAPI_OK)
    std::cerr << "PAPI error " << retval << "\n";
   for (int j = 0; j < num_events_; ++j) {
    counter_bundles_[j].counters.emplace_back(event_counters_[j]);
   }
   PAPI_shutdown();
  }

  void print_header(){
   for (int i = 0; i < num_events_; ++i) {
    PRINT_CSV(counter_bundles_[i].event_name);
   }
  }

  void print(int inst_no){
   for (int i = 0; i < num_events_; ++i) {
    PRINT_CSV(counter_bundles_[i].counters[inst_no]);
   }
  }
 };



}


#endif //FUSION_PAPIWRAPPER_H
