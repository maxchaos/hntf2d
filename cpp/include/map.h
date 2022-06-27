#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>

class HarmonicMap2D
{
private:
  std::vector< Eigen::MatrixXf* > boundaries;
  std::vector< Eigen::VectorXf* > ku, kv;
  Eigen::VectorXf _up, _vp;
public:
  HarmonicMap2D();
  ~HarmonicMap2D();
  Eigen::Matrix2f jacob(float x, float y) const;
  Eigen::Vector2f map(float x, float y) const;
  void print_boundaries() const;
  void boundary_append(const Eigen::MatrixXf &mtx);
  void boundary_set(int index, const Eigen::MatrixXf &mtx);
  Eigen::MatrixXf boundary_get(int index) const;
  void ku_set(int index, const Eigen::VectorXf &vec);
  void kv_set(int index, const Eigen::VectorXf &vec);
  Eigen::VectorXf ku_get(int index) const;
  Eigen::VectorXf kv_get(int index) const;
  void solve(Eigen::VectorXf uo, Eigen::VectorXf vo);
  Eigen::MatrixXf puncts() const;
  size_t boundaries_amount() const { return boundaries.size(); }
  size_t nodes_amount() const;
private:
  std::vector< Eigen::MatrixXf* > _centers;
  std::vector< Eigen::VectorXf* > _lengths;
  std::vector< Eigen::MatrixXf* > _normals;
  bool _state_needs_update;
  void _invalidate_state();
  void _state_update();
private:
  Eigen::MatrixXf _mtx_geom, _mtx_univ, _mtx_aug;
  bool _geom_needs_update, _univ_needs_update, _aug_needs_update;
  void _mtx_geom_build();
  void _mtx_univ_build();
  void _mtx_aug_build();
public:
  Eigen::MatrixXf mtx_geom();
  Eigen::MatrixXf mtx_univ();
  Eigen::MatrixXf mtx_aug();
  /* Eigen::Vector2f boundary(int i, int j); */
  /* void boundary(int i, int j, Eigen::Vector2f &mtx); */
  /* std::vector< MatrixXf* >::iterator boundaries(); */
public:
  void load_xml(const std::string &fn);
  void save_xml(const std::string &fn) const;
};
