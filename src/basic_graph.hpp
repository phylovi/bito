// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <cstddef>
#include <limits>

using VertexId = size_t;
using LineId = size_t;
static constexpr size_t NoId = std::numeric_limits<size_t>::max();

template <typename Traits, template<typename> typename ... Features>
class BasicLine : public Traits, public Features<Traits> ... {
public:

  BasicLine() = default;

  constexpr BasicLine(VertexId first, VertexId second) noexcept : first_{first}, second_{second} {}
    
  constexpr VertexId GetFirst() const noexcept { return first_; }
  constexpr VertexId GetSecond() const noexcept { return second_; }

private:
  template <typename, template<typename> typename ... > friend class BasicGraph;
  
  void OnLineChanged(const BasicLine& line, LineId id) {
    first_ = line.first_;
    second_ = line.second_;
    (Features<Traits>::OnLineChanged(line, id), ...);
  }
  
  VertexId first_ = NoId;
  VertexId second_ = NoId;
};

template <typename Traits, template<typename> typename ... Features>
class BasicVertex : public Traits, public Features<Traits> ... {
public:

  using typename Traits::Line;
  using typename Traits::Neighbors;

  BasicVertex() = default;
  
  constexpr const auto& GetNeighbors() const noexcept { return neighbors_; }

private:
  template <typename, template<typename> typename ... > friend class BasicGraph;
  
  void OnLineChanged(VertexId id, const Line& line, LineId lineId) {
    neighbors_.OnLineChanged(*this, id, line, lineId);
    (Features<Traits>::OnLineChanged(id, line, lineId), ...);
  }

  Neighbors neighbors_;
};

template <typename Traits, template<typename> typename ... Features>
class BasicGraph : public Traits, public Features<Traits> ... {
public:

  using typename Traits::Line;
  using typename Traits::Vertex;
  using typename Traits::Lines;
  using typename Traits::Vertices;

  constexpr const auto& GetLines() const noexcept { return lines_; }
  constexpr const auto& GetVertices() const noexcept { return vertices_; }
  
  const Line& AddLine(LineId id, const Line& newLine) {
    auto& line = GetOrInsert(lines_, id);
    line = newLine;
    OnLineChanged(line, id);
    return line;
  }

  const Vertex& AddVertex(VertexId id, const Vertex& newVertex) {
    auto& vertex = GetOrInsert(vertices_, id);
    vertex = newVertex;
    return vertex;
  }
  
private:
  
  void OnLineChanged(const Line& line, LineId id) {
    GetOrInsert(lines_, id).OnLineChanged(line, id);
    GetOrInsert(vertices_, line.GetFirst()).OnLineChanged(line.GetFirst(), line, id);
    GetOrInsert(vertices_, line.GetSecond()).OnLineChanged(line.GetSecond(), line, id);
  }

  Lines lines_;
  Vertices vertices_;
};
