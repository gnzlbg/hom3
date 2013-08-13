// #ifndef HOM3_FILE_HPP_
// #define HOM3_FILE_HPP_
// ////////////////////////////////////////////////////////////////////////////////
// /// \file \brief This file defines a File handler class for HPC IO
// ////////////////////////////////////////////////////////////////////////////////
// /// Includes:
// #include <type_traits>
// #include <boost/filesystem.hpp>
// #include "../globals.hpp"
// #include <pnetcdf.h>
// /// Options:
// #define ENABLE_DBG_ 0
// #include "../misc/dbg.hpp"
// ////////////////////////////////////////////////////////////////////////////////
// namespace io {
// ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// /// \brief Non-member file handling functions
// namespace file {
// ////////////////////////////////////////////////////////////////////////////////

// /// \brief File/Data type information
// ///@{
// namespace type {
// /// \brief HDF5 file-type tag
// struct HDF5 {};
// /// \brief pnetcdf file-type tag
// struct PNetCDF {};
// }
// /// \brief Output type information
// enum class data_type { INT, FLOAT, DOUBLE, STRING, UINT, LLONG, ULLONG };
// ///@}

// /// \brief IO operation tags
// namespace op {
// /// \brief Read-only operation tag
// struct read {};
// /// \brief Read/Write operation tag
// struct write {};
// }

// /// \brief File information classes
// namespace info {
// /// \brief Type identifiers
// template<class T> using LookUpTable = std::unordered_map<std::string,T>;
// /// \brief Attribute information
// struct Attribute {
//   std::string name;
//   int id;
//   data_type type;
//   Ind length;
// };
// /// \brief Array data information
// struct Array {
//   std::string name;
//   int id;
//   data_type type;
//   std::vector<Ind> dimensions;
//   std::vector<Attribute> attributes;
// };

// }

// /// \brief Does file "fileName" exist?
// bool exists(const std::string fileName, const boost::mpi::communicator& comm) {
//   using namespace boost;
//   return mpi::root_do(comm,[&](){return filesystem::exists(fileName);});
// }
// /// \brief Removes a file _if_ it exists.
// /// \returns True if file existed, false otherwise.
// bool remove(const std::string fileName, const boost::mpi::communicator& comm) {
//   using namespace boost;
//   return mpi::root_do(comm,[&](){return filesystem::remove(fileName);});
// }

// namespace pnetcdf {

// /// \brief Terminates the program in case of error
// void handle_error(int status, const std::string position) {
//   if (status != NC_NOERR) {
//     TERMINATE("ERROR at " + position + "\n PNetCDF error msg: "
//                      + ncmpi_strerror(status));
//   }
// }
// /// \brief Calls a PNetCDF function
// template<class F> void execute(F&& f, const std::string position) {
//   handle_error(f(),position);
// }

// }

// ////////////////////////////////////////////////////////////////////////////////
// } // namespace file
// ////////////////////////////////////////////////////////////////////////////////

// /// \brief File Handler type
// template<class FileType> struct File {};

// /// \brief HDF5 File Handler
// template<> struct File<file::type::HDF5> {};

// /// \brief PNetCDF File Handler
// template<> struct File<file::type::PNetCDF> {

//   //////////////////////////////////////////////////////////////////////////////
//   /// \name Opening / closing files
//   ///@{

//   /// \brief Creates a read-only file handler
//   ///
//   /// Open in read mode:
//   /// - scans datafile for names: attributes, arrays, ...
//   /// - leaves the file open
//   File(const std::string fName, file::op::read,
//        const boost::mpi::communicator& comm = boost::mpi::communicator())
//       : fileName(fName), mpiComm(comm), fileAccessMode(access_mode::read),
//         isOpen_(false) {
//     DBG_ON("reading file");
//     open(file::op::read());
//     build_file_information();
//   }

//   /// \brief Creates a writable file handler
//   File(const std::string fName, file::op::write,
//        const boost::mpi::communicator& comm = boost::mpi::communicator())
//       : fileName (fName), mpiComm(comm), fileAccessMode(access_mode::write),
//         isOpen_(false){}

//   /// \brief Destroyes a PNetCDF file handler (writes all output to disk!)
//   ~File() { close(); };
//   /// \brief File handlers are not default constructible
//   File() = delete;
//   ///@}
//   //////////////////////////////////////////////////////////////////////////////

//   //////////////////////////////////////////////////////////////////////////////
//   /// \name Write/Read file data
//   ///@{
//   /// \brief Returns information about all attributes
//   file::info::LookUpTable<file::info::Attribute> global_attributes() const
//   { return globalAttributesInfo; };
//   /// \brief Reads/Writes attributes
//   template<class T, class M = file::op::read> void attribute
//   (const std::string varName, const std::string attName, T& value, M = M());
//   /// \brief Returns information about all arrays
//   const file::info::LookUpTable<file::info::Array>& arrays() const { return variablesInfo; };
//   /// \brief Reads/Writes array data
//   ///@{
//   template<class T> void array(const std::string n, T *const);
//   //template<class T> void array(const Streamable&);
//   ///@}
//   ///@}
//   //////////////////////////////////////////////////////////////////////////////

//  private:
//   /// File name
//   std::string fileName;
//   /// MPI Communicator (group of MPI Processes that write to file)
//   boost::mpi::communicator mpiComm;

//   /// File access mode
//   ///@{
//   enum class access_mode { read, write };
//   access_mode fileAccessMode;
//   ///@}

//   /// Is a file open or closed?
//   ///@{
//   bool& is_open() { return isOpen_; }
//   bool isOpen_;
//   ///@}

//   /// File identifier
//   int fileId;

//   /// \brief Open/close files
//   ///@{
//   void open(file::op::read) {
//     ASSERT(!is_open(), "File is already open!");
//     DBG_ON("opening file");
//     file::pnetcdf::execute([&](){return ncmpi_open(mpiComm,fileName.c_str(),
//                                   NC_NOWRITE,mpi_info(),&fileId);}, AT_);

//     is_open() = false;
//   }
//   void open(file::op::write) {
//     file::pnetcdf::execute([&](){return ncmpi_create(mpiComm, fileName.c_str(),
//                                     NC_NOCLOBBER|NC_64BIT_OFFSET|NC_64BIT_DATA,
//                                     mpi_info(), &fileId); }, AT_);
//   }
//   void close() {
//     ASSERT(is_open(), "File is not open!");
//     file::pnetcdf::execute([&](){return ncmpi_close(fileId);}, AT_);
//     is_open() = false;
//   }
//   ///@}

//   static const Ind maxAttributeNameLength = 256;

//   /// \name File information
//   ///@{
//   /// \brief Builds file-information look-up tables
//   void build_file_information() {
//     assert_open();

//     int noVariables;
//     file::pnetcdf::execute([&](){ return ncmpi_inq_nvars(fileId, &noVariables); }, AT_);
//     for(int varId = 0; varId < noVariables; ++varId) {
//       char tmpName[NC_MAX_NAME+1];
//       nc_type tmpType;
//       int noDimensions, tmpDimensions[NC_MAX_VAR_DIMS], noAttributes;
//       file::pnetcdf::execute([&](){ return ncmpi_inq_var(fileId,varId,tmpName,&tmpType,
//                                         &noDimensions,tmpDimensions,&noAttributes);
//         }, AT_);

//       std::string name(tmpName);
//       file::data_type type = translate_type(tmpType);
//       std::vector<Ind> dimensions;
//       for(int i = 0; i < noDimensions; ++i) {
//         dimensions.push_back(static_cast<Ind>(tmpDimensions[i]));
//       }
//       std::vector<file::info::Attribute> attributes_;
//       for(int i = 0; i < noAttributes; ++i) {
//         attributes_.push_back(attribute_info(varId,i));
//       }
//       file::info::Array arrayInfo_= {name,varId,type,dimensions,attributes_};
//       variablesInfo.emplace(std::make_pair(name,arrayInfo_));
//     }

//     // Build global attributes lookup table
//     int noGlobalAttributes;
//     file::pnetcdf::execute([&](){ return ncmpi_inq_natts(fileId, &noGlobalAttributes); }, AT_);
//     for(int gattId = 0; gattId < noGlobalAttributes; ++gattId) {
//       auto attribute = attribute_info(NC_GLOBAL,gattId);
//       globalAttributesInfo.emplace(std::make_pair(attribute.name,attribute));
//     }
//   }
//   /// Look-up table of global attributes information
//   file::info::LookUpTable<file::info::Attribute> globalAttributesInfo;
//   /// Look-up table of array variables information
//   file::info::LookUpTable<file::info::Array> variablesInfo;

//   /// \brief Returns attribute info for varId and attId
//   file::info::Attribute attribute_info(int varId, int attId) {
//     assert_open();
//     char tmpName[NC_MAX_NAME+1];
//     nc_type tmpType;
//     MPI_Offset tmpLength;
//     file::pnetcdf::execute([&](){return ncmpi_inq_attname(fileId,varId,attId,tmpName);     }, AT_);
//     file::pnetcdf::execute([&](){return ncmpi_inq_atttype(fileId,varId,tmpName,&tmpType);  }, AT_);
//     file::pnetcdf::execute([&](){return ncmpi_inq_attlen (fileId,varId,tmpName,&tmpLength);}, AT_);
//     return { std::string(tmpName), attId,
//           translate_type(tmpType), static_cast<Ind>(tmpLength) };
//   }
//   /// \brief Translates from NC_TYPE to type::data_type
//   file::data_type translate_type(nc_type ncType) {
//     switch(ncType) {
//       case NC_CHAR:   { return file::data_type::STRING; }
//       case NC_INT:    { return file::data_type::INT;    }
//       case NC_FLOAT:  { return file::data_type::FLOAT;  }
//       case NC_DOUBLE: { return file::data_type::DOUBLE; }
//       case NC_UINT:   { return file::data_type::UINT;   }
//       case NC_INT64:  { return file::data_type::LLONG; }
//       case NC_UINT64: { return file::data_type::ULLONG; }

//       default: { TERMINATE("Unknown nc_type"); }
//     }
//   }
//   ///@}

//   template<class F, class T>
//   void read_attribute
//   (F&& f, const std::string varName, const std::string attName,
//    T&& value, const std::string at_) {
//     if(variablesInfo.find(varName) != end(variablesInfo)) {
//       assert_open();
//       const int varId = variablesInfo[varName].id;
//       auto it = boost::find_if(variablesInfo[varName].attributes,
//                                [&](const file::info::Attribute& a) {
//                                  return a.name == attName;});
//       if(it == std::end(variablesInfo[varName].attributes)) {
//         error::terminate("Variable with name " + varName +
//                          " has no attributed named " + attName, at_);
//       }
//       file::pnetcdf::execute([&](){ return f(fileId,varId,attName.c_str(),&value); }, at_);
//     } else {
//       error::terminate("Variable with name " + varName + " not found", at_);
//     }
//   }


//   void assert_open()   { ASSERT(is_open(),"File is closed!"); }
//   void assert_closed() { ASSERT(!is_open(),"File is open!"); }
// };

// // ////////////////////////////////////////////////////////////////////////////////
// /// \name Attribute functions
// ///@{

// /// \name Readers
// ///@{

// template<> void File<file::type::PNetCDF>::attribute<Ind,file::op::read>
// (const std::string varName, const std::string attName, Ind& value, file::op::read) {
//   static_assert(std::is_same<Ind,unsigned long long>::value, "Needs update!");
//   read_attribute(ncmpi_get_att_ulonglong, varName, attName, value, AT_);
// }

// template<> void File<file::type::PNetCDF>::attribute<Int,file::op::read>
// (const std::string varName, const std::string attName, Int& value, file::op::read) {
//   static_assert(std::is_same<Int,long long>::value, "Needs update!");
//   read_attribute(ncmpi_get_att_longlong, varName, attName, value, AT_);
// }

// template<> void File<file::type::PNetCDF>::attribute<Num,file::op::read>
// (const std::string varName, const std::string attName, Num& value, file::op::read) {
//   static_assert(std::is_same<Num,double>::value, "Needs update!");
//   read_attribute(ncmpi_get_att_double, varName, attName, value, AT_);
// }

// template<> void File<file::type::PNetCDF>::attribute<std::string,file::op::read>
// (const std::string varName, const std::string attName, std::string& value, file::op::read) {
//   char* buf;
//   auto f = [&](int fId, int varId, const char* name, char** val) -> int {
//     MPI_Offset len;
//     file::pnetcdf::execute([&](){ return ncmpi_inq_attlen(fId,varId,name,&len); }, AT_);
//     *val = new char[len+1];
//     return ncmpi_get_att_text(fileId,varId,name,*val);
//   };
//   read_attribute(f, varName, attName, buf, AT_);
//   value = std::string(buf);
//   delete buf;
// }
// ///@}

// /// \name Writers
// ///@{

// // template<> void File<file::type::PNetCDF>::attribute<Ind,op::write>
// // (const std::string name, Ind& value, const std::string dataSet, op::write) {

// // }

// // template<> void File<file::type::PNetCDF>::attribute<Int,op::write>
// // (const std::string name, Int& value, const std::string dataSet, op::write) {

// // }

// // template<> void File<file::type::PNetCDF>::attribute<Num,op::write>
// // (const std::string name, Num& value, const std::string dataSet, op::write) {

// // }

// // template<> void File<file::type::PNetCDF>::attribute<std::string,op::write>
// // (const std::string name, std::string& value, const std::string dataSet, op::write) {

// // }
// ///@}

// // ///@}
// // ////////////////////////////////////////////////////////////////////////////////

// ////////////////////////////////////////////////////////////////////////////////
// } // namespace io
// ////////////////////////////////////////////////////////////////////////////////
// #undef ENABLE_DBG_
// ////////////////////////////////////////////////////////////////////////////////
// #endif
