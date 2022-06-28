#include <hntf2d/map.h>

#include <cmath>
#include <exception>

#include <tinyxml2.h>

HarmonicMap2D::HarmonicMap2D() :
  boundaries{}, ku{}, kv{},
  _state_needs_update{true},
  _geom_needs_update{true},
  _univ_needs_update{true},
  _aug_needs_update{true}
{}

HarmonicMap2D::~HarmonicMap2D()
{
  for( auto elt : this->boundaries )
    if(elt != nullptr)
      delete elt;
  for( auto elt : this->ku )
    if(elt != nullptr)
      delete elt;
  for( auto elt : this->kv )
    if(elt != nullptr)
      delete elt;
}

void
HarmonicMap2D::boundary_append(const Eigen::MatrixXf &geom_)
{
  Eigen::MatrixXf *geom = new Eigen::MatrixXf();
  Eigen::VectorXf *u, *v;
  *geom = geom_;
  u = new Eigen::VectorXf(geom->rows());
  v = new Eigen::VectorXf(geom->rows());
  u->fill(0);
  v->fill(0);
  this->boundaries.push_back(geom);
  this->ku.push_back(u);
  this->kv.push_back(v);
  //// State needs to be updated.
  _invalidate_state();
}

Eigen::MatrixXf
HarmonicMap2D::boundary_get(int index)
  const
{
  return *(this->boundaries.at(index));
}

void
HarmonicMap2D::boundary_set(int index, const Eigen::MatrixXf &mtx)
{
  // if (!(index > -1 && index < this->boundaries.size()))
  //   throw std::out_of_range("boundary index out of range");
  if (mtx.cols() != 2)
    throw std::length_error("invalid matrix shape");
  size_t len_old, len_new;
  len_old = this->boundaries.at(index)->rows();
  len_new = mtx.rows();
  delete (this->boundaries.at(index));
  Eigen::MatrixXf *_mtx = new Eigen::MatrixXf();
  *_mtx = mtx;
  this->boundaries[index] = _mtx;
  if (len_old != len_new)
    {
      delete (this->ku.at(index));
      delete (this->kv.at(index));
      Eigen::VectorXf *_ku, *_kv;
      _ku = new Eigen::VectorXf(_mtx->rows());
      _kv = new Eigen::VectorXf(_mtx->rows());
      _ku->fill(0);
      _kv->fill(0);
      this->ku[index] = _ku;
      this->kv[index] = _kv;
    }
  //// State needs to be updated.
  _invalidate_state();
}

Eigen::VectorXf
HarmonicMap2D::ku_get(int index)
  const
{
  return *(this->ku.at(index));
}

void
HarmonicMap2D::ku_set(int index, const Eigen::VectorXf &vec)
{
  if (vec.size() != ku.at(index)->size())
    throw std::length_error("invalid argument size");
  Eigen::VectorXf *_vec = new Eigen::VectorXf();
  delete (this->ku[index]);
  *_vec = vec;
  this->ku[index] = _vec;
  //// State needs to be updated.
  _invalidate_state();
}

Eigen::VectorXf
HarmonicMap2D::kv_get(int index)
  const
{
  return *(this->kv.at(index));
}

void
HarmonicMap2D::kv_set(int index, const Eigen::VectorXf &vec)
{
  if (vec.size() != kv.at(index)->size())
    throw std::length_error("invalid argument size");
  Eigen::VectorXf *_vec = new Eigen::VectorXf();
  delete (this->kv[index]);
  *_vec = vec;
  this->kv[index] = _vec;
  //// State needs to be updated.
  _invalidate_state();
}

void HarmonicMap2D::print_boundaries()
  const
{
  for (auto elt : this->boundaries)
    std::cout << *elt << std::endl;
}

Eigen::Matrix2f
HarmonicMap2D::jacob( float x, float y )
  const
{
  float xa, ya, xb, yb, xv, yv, xn, yn;
  float xc, yc;
  float x_, y_;
  float xo, yo, xop, xom, xop2, xom2, yo2;
  float l;
  float R[2][2];
  float gradl[2], gradu[2] = {0.0, 0.0}, gradv[2] = {0.0, 0.0};
  Eigen::Matrix2f jacob;
  for( std::size_t i = 0; i < this->boundaries.size(); i++ )
    {
      auto boundary = this->boundaries[i];
      auto &bnd = *boundary;
      auto &ku = *(this->ku[i]);
      auto &kv = *(this->kv[i]);
      for( std::size_t j = 0; j < boundary->rows() - 1; j++ )
        {
          xa = bnd(j,0);
          ya = bnd(j,1);
          xb = bnd(j+1,0);
          yb = bnd(j+1,1);
          xv = xb-xa;
          yv = yb-ya;
          // Length of element.
          l = sqrt(xv*xv+yv*yv);
          xn = xv / l;
          yn = yv / l;
          // Rotation matrix from local to global coordinate frame.
          R[0][0] = xn;
          R[0][1] = -yn;
          R[1][0] = +yn;
          R[1][1] = xn;
          // Center of linear segment.
          xc = 0.5*(xa+xb);
          yc = 0.5*(ya+yb);
          // Point (x,y) w.r.t. to element's local coordinate frame.
          x_ = x-xc;
          y_ = y-yc;
          xo = R[0][0]*x_ + R[1][0]*y_;
          yo = R[0][1]*x_ + R[1][1]*y_;
          // Intermediate values.
          xop = xo+l/2;
          xom = xo-l/2;
          xop2 = xop*xop;
          xom2 = xom*xom;
          yo2 = yo*yo;
          // Unweighted normal gradient expressed in local coordinate frame.
          gradl[0] = -0.5*std::log((yo2 + xom2)/(yo2 + xop2));
          gradl[1] = -std::atan(yo*(xom-xop)/(yo2 + xom*xop));
          // Append element's ij components to gradient of u w.r.t. (x, y).
          gradu[0] += ku(j) * ( R[0][0]*gradl[0] + R[0][1]*gradl[1] );
          gradu[1] += ku(j) * ( R[1][0]*gradl[0] + R[1][1]*gradl[1] );
          // Append element's ij components to gradient of v w.r.t. (x, y).
          gradv[0] += kv(j) * ( R[0][0]*gradl[0] + R[0][1]*gradl[1] );
          gradv[1] += kv(j) * ( R[1][0]*gradl[0] + R[1][1]*gradl[1] );
        }
    }
  jacob(0,0) = gradu[0];
  jacob(0,1) = gradu[1];
  jacob(1,0) = gradv[0];
  jacob(1,1) = gradv[1];
  return jacob;
}

Eigen::Vector2f
HarmonicMap2D::map( float x, float y )
  const
{
  float xa, ya, xb, yb, a, b;
  float axby, aybx, axby2, x2, y2, x2y2, a2b2, xat2, ybt2;
  float xo = x, yo = y;
  float wc;
  float l;
  float u = 0.0, v = 0.0;
  for( std::size_t i = 0; i < this->boundaries.size(); i++ )
    {
      auto boundary = this->boundaries[i];
      auto &bnd = *boundary;
      auto &ku = *(this->ku[i]);
      auto &kv = *(this->kv[i]);
      for( std::size_t j = 0; j < boundary->rows() - 1; j++ )
        {
          xa = bnd(j,0);
          ya = bnd(j,1);
          xb = bnd(j+1,0);
          yb = bnd(j+1,1);
          a = xb-xa;
          b = yb-ya;
          x = xo-xa;
          y = yo-ya;
          axby = a*x+b*y;
          aybx = a*y-b*x;
          axby2 = axby*axby;
          x2y2 = x*x+y*y;
          a2b2 = a*a+b*b;
          // Length of element.
          l = sqrt(a2b2);
          xat2 = x-a; xat2 = xat2*xat2;
          ybt2 = y-b; ybt2 = ybt2*ybt2;
          wc =
            std::log(xat2 + ybt2) +
            axby/a2b2 * std::log(x2y2/(xat2+ybt2)) +
            // (1.0 -axby/a2b2) * std::log(xat2 + ybt2) +
            // axby/a2b2*std::log(x2y2) +
            2.0 * aybx/a2b2 *
            (std::atan(axby/aybx) + std::atan((a2b2-axby)/aybx))
            -2.0;
          wc = l * wc;
          u += ku(j) * wc;
          v += kv(j) * wc;
        }
    }
  return Eigen::Vector2f(u, v);
}

void
HarmonicMap2D::_invalidate_state()
{
  this->_state_needs_update = true;
  this->_geom_needs_update = true;
  this->_univ_needs_update = true;
  this->_aug_needs_update = true;
}

void
HarmonicMap2D::_state_update()
{
  if(!(this->_state_needs_update))
    return;
  for(auto elt : this->_centers) if(elt != nullptr) delete elt;
  for(auto elt : this->_normals) if(elt != nullptr) delete elt;
  for(auto elt : this->_lengths) if(elt != nullptr) delete elt;
  this->_centers.resize(0);
  this->_normals.resize(0);
  this->_lengths.resize(0);
  for(size_t k = 0 ; k < this->boundaries.size(); k++)
    {
      Eigen::MatrixXf *boundary = this->boundaries[k];
      size_t n = boundary->rows() - 1 ;
      Eigen::MatrixXf *centers, *normals;
      Eigen::VectorXf *lengths;
      centers = new Eigen::MatrixXf(n, 2);
      normals = new Eigen::MatrixXf(n, 2);
      lengths = new Eigen::VectorXf(n);
      auto mtx_a = boundary->block(0,0,n,2), mtx_b = boundary->block(1,0,n,2);
      // *centers = 0.5*(boundary->block(1,0,n,2) + boundary->block(0,0,n,2));
      *centers = 0.5*(mtx_a + mtx_b);
      // Eigen::MatrixXf delta =
      //        boundary->block(1,0,n,2) - boundary->block(0,0,n,2);
      // *lengths = (delta * delta).rowwise().sum();
      // *lengths = (boundary->block(1,0,n,2) - boundary->block(0,0,n,2)).rowwise().squaredNorm();
      Eigen::MatrixXf deltas(mtx_b - mtx_a);
      if(k == 0)
        {
          normals->block(0,0,n,1) = -deltas.block(0,1,n,1);
          normals->block(0,1,n,1) = deltas.block(0,0,n,1);
        }
      else
        {
          normals->block(0,0,n,1) = deltas.block(0,1,n,1);
          normals->block(0,1,n,1) = -deltas.block(0,0,n,1);
        }
      // *lengths = (mtx_b - mtx_a).rowwise().norm();
      *lengths = (deltas).rowwise().norm();
      //// Normalize normals...
      for(int l = 0; l < normals->rows(); l++) {
        double len = (*lengths)(l);
        (*normals)(l,0) /= len;
        (*normals)(l,1) /= len;
      }
      // this->_centers[k] = new Eigen::MatrixXf(n, 2);
      // this->_normals[k] = new Eigen::MatrixXf(n, 2);
      // this->_lengths[k] = new Eigen::VectorXf(n);
      this->_centers.push_back(centers);
      this->_normals.push_back(normals);
      this->_lengths.push_back(lengths);
    }
  this->_state_needs_update = false;
}

void
HarmonicMap2D::_mtx_geom_build()
{
  if(this->_state_needs_update)
    _state_update();
  //// Totals number of elements.
  size_t n = 0;
  for(auto elt : this->_centers)
    n += elt->rows();
  size_t ncs[this->boundaries.size() + 1];
  ncs[0] = 0;
  for(int k = 0; k < this->_centers.size(); k++)
    ncs[k+1] = ncs[k] + this->_centers[k]->rows();
  //// Allocate geometry matrix and fill with zeros.
  this->_mtx_geom.resize(n,n);
  this->_mtx_geom.fill(0);

  for(int i = 0; i < this->boundaries.size(); i++)
    {
      Eigen::MatrixXf *boundary = this->boundaries[i];
      for(int j = 0; j < boundary->rows()-1; j++)
        {
          float xa = (*boundary)(j,0), ya = (*boundary)(j,1);
          float xb = (*boundary)(j+1,0), yb = (*boundary)(j+1,1);
          float a = xb - xa, b = yb - ya;
          float a2b2 = a*a + b*b;
          float len = std::sqrt(a2b2);
          for(int k = 0; k < this->_centers.size(); k++)
            {
              Eigen::MatrixXf *points = this->_centers[k];
              Eigen::VectorXf xo = points->col(0);
              Eigen::VectorXf yo = points->col(1);
              Eigen::VectorXf wc(points->rows());
              for(int l = 0; l < points->rows(); l++)
                {
                  float x = xo(l) - xa, y = yo(l) - ya;
                  float axby = a*x + b*y, aybx = std::abs(a*y - b*x);
                  // float axby2 = axby*axby;
                  float x2 = x*x, y2 = y*y;
                  float x2y2 = x2 + y2;
                  float xat2 = (x-a)*(x-a), ybt2 = (y-b)*(y-b);
                  wc(l) =
                    std::log(xat2 + ybt2) +
                    axby/a2b2 * std::log(x2y2 / (xat2+ybt2)) +
                    2.0 * aybx/a2b2 *
                    (std::atan2(axby, aybx) + std::atan2(a2b2-axby,aybx))
                    -2.0;
                  wc(l) *= len;
                }
              this->_mtx_geom.block(ncs[k], ncs[i]+j, ncs[k+1]-ncs[k], 1) = wc;
            }
        }
    }
}

Eigen::MatrixXf
HarmonicMap2D::mtx_geom()
{
  if(this->_geom_needs_update)
    _mtx_geom_build();
  return this->_mtx_geom;
}

void
HarmonicMap2D::_mtx_univ_build()
{
  if(this->_state_needs_update)
    _state_update();
  //// Totals number of elements.
  size_t n = 0;
  for(auto elt : this->_centers)
    n += elt->rows();
  size_t ncs[this->boundaries.size() + 1];
  ncs[0] = 0;
  for(int k = 0; k < this->_centers.size(); k++)
    ncs[k+1] = ncs[k] + this->_centers[k]->rows();
  //// Allocate geometry matrix and fill with zeros.
  this->_mtx_univ.resize(this->boundaries.size() - 1,n);
  this->_mtx_univ.fill(0);
  if (this->_mtx_univ.rows() == 0)
    return;
  for(int i = 0; i < this->boundaries.size(); i++)
    {
      Eigen::MatrixXf *boundary = this->boundaries[i];
      for(int j = 0; j < boundary->rows()-1; j++)
        {
          float xa = (*boundary)(j,0), ya = (*boundary)(j,1);
          float xb = (*boundary)(j+1,0), yb = (*boundary)(j+1,1);
          float a = xb - xa, b = yb - ya;
          float a2b2 = a*a + b*b;
          float len = std::sqrt(a2b2);
          Eigen::Matrix2f R;
          R << xb-xa, -yb+ya, yb-ya, xb-xa;
          R /= len;
          Eigen::Vector2f pc = this->_centers[i]->row(j);
          for(int k = 1; k < this->_centers.size(); k++)
            {
              Eigen::MatrixXf *points = this->_centers[k];
              // Eigen::MatrixXf po = (points->rowwise() - pc.transpose()) * R;
              for(int l = 0; l < points->rows(); l++)
                {
                  Eigen::Vector2f pt = points->row(l).transpose() - pc;
                  Eigen::Vector2f po = R.transpose() * pt;
                  float xo = po(0), yo = po(1);
                  float xop = xo + len/2, xom = xo - len/2;
                  float xop2 = xop*xop, xom2 = xom*xom, yo2 = yo*yo;
                  Eigen::Vector2f gradl(
                    std::log((yo2+xom2)/(yo2+xop2))/2,
                    std::atan(xom/yo) - std::atan(xop/yo)
                  );
                  Eigen::Vector2f gradw = R * gradl;
                  Eigen::Vector2f normal = this->_normals[k]->row(l);
                  float delta =
                    (gradw.transpose() * normal).sum();
                  if(i == k and j == l)
                    delta = std::abs(delta);
                  this->_mtx_univ(k-1,ncs[i]+j) += delta;
                }
            }
        }
    }
}

Eigen::MatrixXf
HarmonicMap2D::mtx_univ()
{
  if(this->_univ_needs_update)
    _mtx_univ_build();
  return this->_mtx_univ;
}

void
HarmonicMap2D::_mtx_aug_build()
{
  if(this->_state_needs_update)
    _state_update();
  if(this->_geom_needs_update)
    _mtx_geom_build();
  if(this->_univ_needs_update)
    _mtx_univ_build();

  size_t ne = this->_mtx_geom.rows();
  size_t np = this->_mtx_univ.rows();
  size_t nt = ne + np;

  if(np == 0)
    this->_mtx_aug = this->_mtx_geom;
  else
    {
      float agm = this->_mtx_geom.array().abs().maxCoeff();
      float aum = this->_mtx_univ.array().abs().maxCoeff();
      this->_mtx_aug.resize(nt, nt);
      this->_mtx_aug.block(0,0,ne,ne) = this->_mtx_geom;
      this->_mtx_aug.block(ne,0,np,ne) = agm/aum * (this->_mtx_univ);
      this->_mtx_aug.block(0,ne,nt,np).fill(0);

      size_t n = 0;
      for(auto elt : this->_centers)
        n += elt->rows();

      size_t ncs[this->boundaries.size() + 1];
      ncs[0] = 0;
      for(int k = 0; k < this->_centers.size(); k++)
        ncs[k+1] = ncs[k] + this->_centers[k]->rows();

      for(int k = 1, idx = this->_centers[0]->rows(), len = 0;
          k < this->boundaries.size();
          k++)
        {
          len = this->_centers[k]->rows();
          this->_mtx_aug.block(idx,ne+k-1,len,1).fill(-1);
          idx += len;
        }
    }
}

Eigen::MatrixXf
HarmonicMap2D::mtx_aug()
{
  if(this->_aug_needs_update)
    _mtx_aug_build();
  return this->_mtx_aug;
}

void
HarmonicMap2D::solve(Eigen::VectorXf uo, Eigen::VectorXf vo)
{
  if(this->_aug_needs_update)
    _mtx_aug_build();
  if (uo.rows() != this->boundaries[0]->rows()-1)
    throw std::domain_error("unexpected length of uo");
  if (vo.rows() != this->boundaries[0]->rows()-1)
    throw std::domain_error("unexpected length of vo");
  Eigen::VectorXf bu(this->_mtx_aug.rows());
  Eigen::VectorXf bv(this->_mtx_aug.rows());
  bu.block(0,0,uo.rows(),1) = uo;
  bv.block(0,0,vo.rows(),1) = vo;
  int pad_len = this->_mtx_aug.rows() - (this->boundaries[0]->rows() - 1);
  bu.block(uo.rows(),0,pad_len,1).fill(0);
  bv.block(vo.rows(),0,pad_len,1).fill(0);
  Eigen::HouseholderQR< Eigen::MatrixXf > solver(this->_mtx_aug);
  Eigen::VectorXf ku = solver.solve(bu);
  Eigen::VectorXf kv = solver.solve(bv);
  int idx = 0;
  for(int k = 0; k < this->boundaries.size(); k++)
    {
      int len = this->boundaries[k]->rows() - 1;
      *(this->ku[k]) = ku.block(idx,0,len,1);
      *(this->kv[k]) = kv.block(idx,0,len,1);
      idx += len;
    }
  this->_up = ku.block(idx,0,this->boundaries.size()-1,1);
  this->_vp = kv.block(idx,0,this->boundaries.size()-1,1);
}

Eigen::MatrixXf
HarmonicMap2D::puncts()
  const
{
  Eigen::MatrixXf res(this->_up.rows(), 2);
  res << this->_up, this->_vp;
  return res;
}

size_t
HarmonicMap2D::nodes_amount()
  const
{
  size_t res = 0;
  for(size_t k = 0; k < this->boundaries.size(); k++ )
    res += boundaries[k]->rows();
  return res;
}

void
HarmonicMap2D::save_xml(const std::string &fn)
  const
{
  tinyxml2::XMLDocument xml_doc;
  auto xml_root = xml_doc.NewElement("map");
  xml_doc.InsertFirstChild(xml_root);
  // xml_root->SetAttribute("genus", int(this->boundaries.size()-1));
  char buff[30];
  auto num2str = [&](float num) { sprintf(buff, "%e", num); return buff; };
  // Insert boundaries in order.
  for (size_t idx = 0; idx < this->boundaries_amount(); idx++)
    {
      // Aliases for fast access.
      const Eigen::MatrixXf &bnd = *(this->boundaries[idx]);
      const Eigen::VectorXf &ku = *(this->ku[idx]);
      const Eigen::VectorXf &kv = *(this->kv[idx]);
      const Eigen::MatrixXf &centers = *(this->_centers[idx]);
      const Eigen::VectorXf &lengths = *(this->_lengths[idx]);
      const Eigen::MatrixXf &normals = *(this->_normals[idx]);
      // Create node for corresponding boundary.
      auto xml_bnd = xml_doc.NewElement("boundary");
      xml_root->InsertEndChild(xml_bnd);
      // TODO
      // Create node for boundary's actual representation.
      auto xml_pgn = xml_doc.NewElement("curve");
      xml_bnd->InsertEndChild(xml_pgn);
      // xml_pgn->SetAttribute("size", bnd.rows());
      // Store representation.
      for (size_t i = 0; i < bnd.rows(); i++)
        {
          auto xml_point = xml_doc.NewElement("point");
          xml_point->SetAttribute("x", num2str(bnd(i, 0)));
          xml_point->SetAttribute("y", num2str(bnd(i, 1)));
          xml_pgn->InsertEndChild(xml_point);
        }
      // Create node for list of boundary's elements.
      auto xml_elts = xml_doc.NewElement("elements");
      xml_bnd->InsertEndChild(xml_elts);
      // xml_elts->SetAttribute("size", centers.rows());
      for (size_t i = 0; i < centers.rows(); i++)
        {
          auto xml_elt = xml_doc.NewElement("element");
          xml_elt->SetAttribute("length", lengths(i));
          xml_elt->SetAttribute("cx", num2str(centers(i,0)));
          xml_elt->SetAttribute("cy", num2str(centers(i,1)));
          xml_elt->SetAttribute("nx", num2str(normals(i,0)));
          xml_elt->SetAttribute("ny", num2str(normals(i,1)));
          xml_elt->SetAttribute("ku", num2str(ku(i)));
          xml_elt->SetAttribute("kv", num2str(kv(i)));
          xml_elts->InsertEndChild(xml_elt);
        }
      // Store associated puncture.
      if(idx > 0) {
        auto xml_punct = xml_doc.NewElement("puncture");
        xml_bnd->InsertEndChild(xml_punct);
        xml_punct->SetAttribute("u", num2str(this->_up(idx-1)));
        xml_punct->SetAttribute("v", num2str(this->_vp(idx-1)));
      }
    }
  // // Create node for list of punctures.
  // // auto xml_puncts = xml_doc.NewElement("punctures");
  // // xml_root->InsertEndChild(xml_puncts);
  // // Store punctures.
  // for(size_t i = 0; i < this->_up.rows(); i++)
  //   {
  //     auto xml_punct = xml_doc.NewElement("puncture");
  //     xml_root->InsertEndChild(xml_punct);
  //     // xml_puncts->InsertEndChild(xml_punct);
  //     xml_punct->SetAttribute("u", this->_up(i));
  //     xml_punct->SetAttribute("v", this->_vp(i));
  //   }
  xml_doc.SaveFile(fn.c_str());
}

void
HarmonicMap2D::load_xml(const std::string &fn)
{
  tinyxml2::XMLDocument xml_doc;
  {
    auto res = xml_doc.LoadFile(fn.c_str());
    if (res != tinyxml2::XML_SUCCESS)
      throw "failed to load harmonic map from xml file";
  }
  auto xml_root = xml_doc.RootElement();
  //// Load xml boundaries.
  std::vector< tinyxml2::XMLElement* > xml_bnds;
  {
    auto xml_bnd = xml_root->FirstChildElement("boundary");
    if(xml_bnd == nullptr)
      throw "map has no 'boundary' elements";
    while (xml_bnd != nullptr) {
      xml_bnds.push_back(xml_bnd);
      xml_bnd = xml_bnd->NextSiblingElement("boundary");
    }
  }
  //// Infer domain's genus.
  const size_t genus = xml_bnds.size() - 1;
  this->_up.resize(genus);
  this->_vp.resize(genus);
  //// Load each individual boundary.
  for(size_t bnd_idx = 0; bnd_idx < xml_bnds.size(); bnd_idx++)
    {
      auto xml_bnd = xml_bnds[bnd_idx];
      // Load boundary's polygonal representation.
      {
        auto xml_pgn = xml_bnd->FirstChildElement("curve");
        if (xml_pgn == nullptr)
          throw "boundary has no 'curve' element";
        std::vector<tinyxml2::XMLElement*> xml_pnts;
        // Load xml points.
        {
          auto xml_pnt = xml_pgn->FirstChildElement("point");
          if (xml_pnt == nullptr)
            throw "curve has no 'point' elements";
          while (xml_pnt != nullptr) {
            xml_pnts.push_back(xml_pnt);
            xml_pnt = xml_pnt->NextSiblingElement("point");
          }
        }
        // Infer polygon's size.
        const size_t pgn_size = xml_pnts.size();
        // Allocate polygon.
        Eigen::MatrixXf *pgn = new Eigen::MatrixXf(pgn_size, 2);
        // Set polygon.
        for (size_t pnt_idx = 0; pnt_idx < pgn_size; pnt_idx++) {
          auto xml_pnt = xml_pnts[pnt_idx];
          //// Read x coordinate.
          if (not xml_pnt->Attribute("x"))
            throw "'point' element has no 'x' attribute";
          (*pgn)(pnt_idx, 0) = xml_pnt->DoubleAttribute("x");
          //// Read y coordinate.
          if (not xml_pnt->Attribute("y"))
            throw "'point' element has no 'y' attribute";
          (*pgn)(pnt_idx, 1) = xml_pnt->DoubleAttribute("y");
        }
        this->boundaries.push_back(pgn);
      }
      //// Load boundary's elements.
      {
        auto xml_elements = xml_bnd->FirstChildElement("elements");
        if (xml_elements == nullptr)
          throw "boundary has no 'elements' element";
        //// Load xml elements.
        std::vector<tinyxml2::XMLElement*> xml_elts;
        {
          auto xml_elt = xml_elements->FirstChildElement("element");
          if (xml_elt == nullptr)
            throw "elements has no 'element' element";
          while(xml_elt != nullptr) {
            xml_elts.push_back(xml_elt);
            xml_elt = xml_elt->NextSiblingElement("element");
          }
        }
        //// Infer amount of elements.
        const size_t elts_amount = xml_elts.size();
        //// Allocate matrices.
        Eigen::MatrixXf *centers = new Eigen::MatrixXf(elts_amount, 2);
        Eigen::MatrixXf *normals = new Eigen::MatrixXf(elts_amount, 2);
        Eigen::VectorXf *lengths = new Eigen::VectorXf(elts_amount);
        Eigen::VectorXf *ku = new Eigen::VectorXf(elts_amount);
        Eigen::VectorXf *kv = new Eigen::VectorXf(elts_amount);
        // while(xml_elt != nullptr)
        for(size_t elt_idx = 0; elt_idx < elts_amount; elt_idx++)
          {
            auto xml_elt = xml_elts[elt_idx];
            //// Read x coordinate of element's center.
            if (not xml_elt->Attribute("cx"))
              throw "'element' has no 'cx' attribute";
            (*centers)(elt_idx, 0) = xml_elt->DoubleAttribute("cx");
            //// Read y coordinate of element's center.
            if (not xml_elt->Attribute("cy"))
              throw "'element' has no 'cy' attribute";
            (*centers)(elt_idx, 1) = xml_elt->DoubleAttribute("cy");
            //// Read x coordinate of element's normal.
            if (not xml_elt->Attribute("nx"))
              throw "'element' has no 'nx' attribute";
            (*normals)(elt_idx, 0) = xml_elt->DoubleAttribute("nx");
            //// Read y coordinate of element's normal.
            if (not xml_elt->Attribute("ny"))
              throw "'element' has no 'ny' attribute";
            (*normals)(elt_idx, 1) = xml_elt->DoubleAttribute("ny");
            //// Read element's length.
            if (not xml_elt->Attribute("length"))
              throw "'element' has no 'length' attribute";
            (*lengths)(elt_idx) = xml_elt->DoubleAttribute("length");
            //// Read element's k_u.
            if (not xml_elt->Attribute("ku"))
              throw "'element' has no 'ku' attribute";
            (*ku)(elt_idx) = xml_elt->DoubleAttribute("ku");
            //// Read element's k_v.
            if (not xml_elt->Attribute("kv"))
              throw "'element' has no 'kv' attribute";
            (*kv)(elt_idx) = xml_elt->DoubleAttribute("kv");
          }
        //// Store element parameters.
        this->_centers.push_back(centers);
        this->_normals.push_back(normals);
        this->_lengths.push_back(lengths);
        this->ku.push_back(ku);
        this->kv.push_back(kv);
      }
      //// Load corresponding puncture.
      {
        if (bnd_idx > 0) {
          auto xml_punct = xml_bnd->FirstChildElement("puncture");
          if (xml_punct == nullptr)
            throw "boundary has no 'puncture' element";
          //// Read puncture's u coordinate.
          if (not xml_punct->Attribute("u"))
            throw "puncture has no 'u' attribute";
          this->_up(bnd_idx-1) = xml_punct->DoubleAttribute("u");
          //// Read puncture's v coordinate.
          if (not xml_punct->Attribute("v"))
            throw "puncture has no 'v' attribute";
          this->_vp(bnd_idx-1) = xml_punct->DoubleAttribute("v");
        }
      }
    }
}
