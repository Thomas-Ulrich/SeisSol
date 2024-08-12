#include "BinaryWriter.hpp"
#include <IO/Writer/Instructions/Binary.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "utils/logger.h"

namespace seissol::io::writer::file {

BinaryFile::BinaryFile(MPI_Comm comm) : comm(comm) {}
void BinaryFile::openFile(const std::string& name) {
  MPI_File_open(comm,
                name.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
                MPI_INFO_NULL,
                &file);
}
void BinaryFile::writeGlobal(const void* data, std::size_t size) {
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    MPI_File_write(file, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(comm);
}
void BinaryFile::writeDistributed(const void* data, std::size_t size) {
  // TODO: get max write size
  MPI_File_write_all(file, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
}
void BinaryFile::closeFile() { MPI_File_close(&file); }

BinaryWriter::BinaryWriter(MPI_Comm comm) : comm(comm) {}

void BinaryWriter::write(const async::ExecInfo& info, const instructions::BinaryWrite& write) {
  if (openFiles.find(write.filename) == openFiles.end()) {
    openFiles[write.filename] = std::make_unique<BinaryFile>(BinaryFile(comm));
    openFiles[write.filename]->openFile(write.filename);
  }

  const void* dataPointer = write.dataSource->getPointer(info);

  // TODO: add dimensions
  const auto dataSize = write.dataSource->count(info) * write.dataSource->datatype()->size();

  if (write.dataSource->distributed()) {
    openFiles[write.filename]->writeDistributed(write.dataSource->getPointer(info), dataSize);
  } else {
    openFiles[write.filename]->writeGlobal(write.dataSource->getPointer(info), dataSize);
  }
}

void BinaryWriter::finalize() {
  for (auto& [_, file] : openFiles) {
    file->closeFile();
  }
}

} // namespace seissol::io::writer::file
