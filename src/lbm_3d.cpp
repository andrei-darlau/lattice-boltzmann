// Single-threaded D3Q19 lattice Boltzmann implementation.
// Flow past a generic stationary ogive-nose / cylindrical / boattail body.
// This geometry is not a reproduction of a projectile or other weapon.
// All quantities are in non-dimensional lattice units.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
constexpr int Q = 19;
constexpr std::array<int, Q> cx = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0};
constexpr std::array<int, Q> cy = {0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1};
constexpr std::array<int, Q> cz = {0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1};
constexpr std::array<int, Q> opposite = {0, 2,1, 4,3, 6,5, 8,7, 10,9, 12,11, 14,13, 16,15, 18,17};
constexpr std::array<double, Q> weight = {
    1.0/3.0,
    1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,
    1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,
    1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

struct Config {
  int nx = 48, ny = 48, nz = 112;
  int steps = 500, every = 250;
  double u_in = 0.045;
  double reynolds = 35.0; // nominal Re using full maximum body diameter
  double smagorinsky = 0.0; // experimental eddy-viscosity option, off by default
  std::string out = "output_3d";
};
struct Macros { double rho, ux, uy, uz; };

void usage(const char* name) {
  std::cout << "Usage: " << name << " [--nx N] [--ny N] [--nz N]"
               " [--steps N] [--every N] [--u U] [--re RE]"
               " [--smagorinsky CS] [--out DIRECTORY]\n"
               "Serial D3Q19 LBM around a GENERIC ogive/cylinder/boattail solid, +Z flow.\n";
}
Config parse_args(int argc, char** argv) {
  Config c;
  for (int k=1; k<argc; ++k) {
    std::string a=argv[k];
    if (a=="--help" || a=="-h") { usage(argv[0]); std::exit(0); }
    if (k+1==argc) throw std::invalid_argument("Missing argument after " + a);
    const std::string v=argv[++k];
    if (a=="--nx") c.nx=std::stoi(v);
    else if (a=="--ny") c.ny=std::stoi(v);
    else if (a=="--nz") c.nz=std::stoi(v);
    else if (a=="--steps") c.steps=std::stoi(v);
    else if (a=="--every") c.every=std::stoi(v);
    else if (a=="--u") c.u_in=std::stod(v);
    else if (a=="--re") c.reynolds=std::stod(v);
    else if (a=="--smagorinsky") c.smagorinsky=std::stod(v);
    else if (a=="--out") c.out=v;
    else throw std::invalid_argument("Unknown option: " + a);
  }
  if (c.nx<28 || c.ny<28 || c.nz<64 || c.steps<1 || c.every<1)
    throw std::invalid_argument("Grid must be at least 28x28x64; steps/every >= 1");
  if (!(c.u_in>0 && c.u_in<0.12) || !(c.reynolds>0) ||
      !(c.smagorinsky>=0.0 && c.smagorinsky<=0.2))
    throw std::invalid_argument("Require 0<u<0.12, Re>0, 0<=Cs<=0.2");
  if (c.out.empty()) throw std::invalid_argument("Output directory must not be empty");
  return c;
}
inline std::size_t idx(int x,int y,int z,int nx,int ny) {
  return (static_cast<std::size_t>(z)*ny + y)*nx + x;
}
inline double feq(int q, double rho, double ux, double uy, double uz) {
  const double cu = cx[q]*ux + cy[q]*uy + cz[q]*uz;
  const double uu = ux*ux + uy*uy + uz*uz;
  return weight[q]*rho*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*uu);
}
inline Macros macros(const double* f) {
  double rho=0, mx=0, my=0, mz=0;
  for (int q=0;q<Q;++q) {
    rho+=f[q]; mx+=cx[q]*f[q]; my+=cy[q]*f[q]; mz+=cz[q]*f[q];
  }
  if (!std::isfinite(rho) || !(rho>0.0))
    throw std::runtime_error("Nonpositive/nonfinite density: run is unstable");
  return {rho,mx/rho,my/rho,mz/rho};
}
// The ogive profile is defined in lattice units, not physical dimensions.
// s = distance from the front tip; 0 <= s <= nose_length.
// The formula gives R(0)=0 and R(L)=radius, with zero slope at R(L).
inline double tangent_ogive(double s, double nose_length, double radius) {
  if (s <= 0.0) return 0.0;
  if (s >= nose_length) return radius;
  const double rho = (nose_length*nose_length + radius*radius)/(2.0*radius);
  return std::sqrt(std::max(0.0, rho*rho - (nose_length-s)*(nose_length-s)))
         + radius - rho;
}

// Generic surface of revolution used directly by the LBM solid-mask builder
struct RevolvedProfile {
  double tip, nose_end, barrel_end, tail_end, max_radius, end_radius;

  double radius_at(double z) const {
    if (z < tip || z > tail_end) return 0.0;
    if (z < nose_end) return tangent_ogive(z-tip, nose_end-tip, max_radius);
    if (z < barrel_end) return max_radius;
    const double t=(z-barrel_end)/(tail_end-barrel_end);
    return max_radius + t*(end_radius-max_radius);
  }
};

class Solver {
 public:
  explicit Solver(Config c) : cfg_(std::move(c)),
      n_(static_cast<std::size_t>(cfg_.nx)*cfg_.ny*cfg_.nz),
      solid_(n_,0), f_(n_*Q,0), post_(n_*Q,0), next_(n_*Q,0) {
    radius_=std::min(cfg_.nx,cfg_.ny)/8.5;
    const double diameter=2.0*radius_;
    viscosity_=cfg_.u_in*diameter/cfg_.reynolds;
    tau_=0.5+3.0*viscosity_;
    if (tau_<0.515)
      throw std::invalid_argument("tau too near 0.5 for this educational BGK solver; lower Re");

    // Entire body is defined by a continuous radial profile around the Z axis:
    //   x = cx + R(z)*cos(phi), y = cy + R(z)*sin(phi).
    // Only the inside/outside test is needed for the Cartesian D3Q19 lattice.
    // These arbitrary dimensionless ratios describe a generic test body.
    const double start=0.16*cfg_.nz;
    const double nose_length=2.8*radius_;
    const double mid_length=2.0*radius_;
    const double tail_length=2.7*radius_;
    profile_={start,
              start+nose_length,
              start+nose_length+mid_length,
              start+nose_length+mid_length+tail_length,
              radius_, 0.60*radius_};
    if (profile_.tail_end>0.77*cfg_.nz)
      throw std::invalid_argument("Shape too long for streamwise grid; increase --nz");
    const double xc=(cfg_.nx-1.0)/2.0, yc=(cfg_.ny-1.0)/2.0;
    std::size_t solid_count=0;
    for(int z=0;z<cfg_.nz;++z) {
      const double r=profile_.radius_at(static_cast<double>(z));
      for(int y=0;y<cfg_.ny;++y) for(int x=0;x<cfg_.nx;++x) {
        const auto i=idx(x,y,z,cfg_.nx,cfg_.ny);
        const double dx=x-xc, dy=y-yc;
        solid_[i]=(r>0.0 && dx*dx+dy*dy<=r*r) ? 1 : 0;
        solid_count+=solid_[i];
        for (int q=0;q<Q;++q)
          f_[i*Q+q]=feq(q,1.0,0.0,0.0,cfg_.u_in);
      }
    }
    if (!solid_count) throw std::runtime_error("No solid cells in geometry");
  }

  void run() {
    std::filesystem::create_directories(cfg_.out);
    std::cout << "Serial D3Q19 BGK, grid " << cfg_.nx << 'x' << cfg_.ny << 'x' << cfg_.nz
              << ", generic ogive body radius=" << radius_
              << ", U=" << cfg_.u_in << ", Re=" << cfg_.reynolds
              << ", nu=" << viscosity_ << ", tau=" << tau_
              << ", Cs=" << cfg_.smagorinsky << "\n";
    write_vtk(0);
    for (int t=1;t<=cfg_.steps;++t) {
      collide();
      stream();
      set_open_boundaries();
      f_.swap(next_);
      if (t%cfg_.every==0 || t==cfg_.steps) {
        report(t);
        write_vtk(t);
      }
    }
  }
 private:
  Config cfg_;
  std::size_t n_;
  std::vector<std::uint8_t> solid_;
  std::vector<double> f_, post_, next_;
  double radius_=0, viscosity_=0, tau_=0;
  RevolvedProfile profile_{};

  void collide() {
    for (std::size_t i=0;i<n_;++i) {
      if (solid_[i]) continue;
      const double* fi=&f_[i*Q];
      const Macros m=macros(fi);
      double local_tau=tau_;
      if (cfg_.smagorinsky>0) {
        // Estimate the non-equilibrium momentum-flux norm. Experimental
        // Smagorinsky closure; not a verified turbulence model.
        double pxx=0,pyy=0,pzz=0,pxy=0,pxz=0,pyz=0;
        for (int q=0;q<Q;++q) {
          const double neq=fi[q]-feq(q,m.rho,m.ux,m.uy,m.uz);
          pxx+=cx[q]*cx[q]*neq; pyy+=cy[q]*cy[q]*neq;
          pzz+=cz[q]*cz[q]*neq; pxy+=cx[q]*cy[q]*neq;
          pxz+=cx[q]*cz[q]*neq; pyz+=cy[q]*cz[q]*neq;
        }
        const double norm=std::sqrt(pxx*pxx+pyy*pyy+pzz*pzz+
                                    2*(pxy*pxy+pxz*pxz+pyz*pyz));
        const double cs2=cfg_.smagorinsky*cfg_.smagorinsky;
        local_tau=0.5*(tau_+std::sqrt(tau_*tau_+
                      36.0/std::sqrt(2.0)*cs2*norm/m.rho));
      }
      const double omega=1.0/local_tau;
      double* dst=&post_[i*Q];
      for (int q=0;q<Q;++q) {
        const double eq=feq(q,m.rho,m.ux,m.uy,m.uz);
        dst[q]=fi[q]-omega*(fi[q]-eq);
      }
    }
  }
  void stream() {
    const int nx=cfg_.nx,ny=cfg_.ny,nz=cfg_.nz;
    // Interior +Z flow planes; the Z inlet/outlet are set below.
    // X/Y transverse boundaries are periodic.
    for (int z=1;z<nz-1;++z) for(int y=0;y<ny;++y)
      for (int x=0;x<nx;++x) {
        const auto dest=idx(x,y,z,nx,ny);
        if (solid_[dest]) continue;
        double* dst=&next_[dest*Q];
        for (int q=0;q<Q;++q) {
          const int sx=(x-cx[q]+nx)%nx;
          const int sy=(y-cy[q]+ny)%ny;
          const int sz=z-cz[q];
          const auto src=idx(sx,sy,sz,nx,ny);
          dst[q]=solid_[src] ? post_[dest*Q+opposite[q]]
                              : post_[src*Q+q];
        }
      }
  }
  void set_open_boundaries() {
    // First-order non-equilibrium extrapolation at Z inlet/outlet.
    // This is a simplified educational boundary condition.
    const int nx=cfg_.nx,ny=cfg_.ny,nz=cfg_.nz;
    for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) {
      const auto in=idx(x,y,0,nx,ny);
      const auto in_adj=idx(x,y,1,nx,ny);
      const Macros a=macros(&next_[in_adj*Q]);
      for (int q=0;q<Q;++q)
        next_[in*Q+q]=feq(q,a.rho,0,0,cfg_.u_in)+
          next_[in_adj*Q+q]-feq(q,a.rho,a.ux,a.uy,a.uz);
      const auto out=idx(x,y,nz-1,nx,ny);
      const auto out_adj=idx(x,y,nz-2,nx,ny);
      const Macros b=macros(&next_[out_adj*Q]);
      for (int q=0;q<Q;++q)
        next_[out*Q+q]=feq(q,1.0,b.ux,b.uy,b.uz)+
          next_[out_adj*Q+q]-feq(q,b.rho,b.ux,b.uy,b.uz);
    }
  }

  void report(int t) const {
    double minrho=std::numeric_limits<double>::infinity();
    double maxrho=0, maxu=0;
    for (std::size_t i=0;i<n_;++i) {
      if (solid_[i]) continue;
      const auto m=macros(&f_[i*Q]);
      minrho=std::min(minrho,m.rho);
      maxrho=std::max(maxrho,m.rho);
      maxu=std::max(maxu,std::sqrt(m.ux*m.ux+m.uy*m.uy+m.uz*m.uz));
    }
    if (minrho<0.3 || maxrho>2.0 || maxu>0.35)
      throw std::runtime_error("Flow became unstable; reduce U/Re or refine grid");
    std::cout << "step " << std::setw(6) << t << " | rho ["
              << std::fixed << std::setprecision(5)
              << minrho << ", " << maxrho << "] | max |u|=" << maxu << '\n';
  }

  // VTK XML appended data uses explicitly little-endian byte counts and values.
  static void little_endian(std::ostream& out, std::uint64_t value, int bytes) {
    for (int i=0; i<bytes; ++i)
      out.put(static_cast<char>((value >> (8*i)) & 0xff));
  }
  static void write_float(std::ostream& out, float value) {
    static_assert(sizeof(float)==4 && std::numeric_limits<float>::is_iec559,
                  "VTK Float32 requires IEEE 754 float");
    std::uint32_t bits;
    std::memcpy(&bits,&value,sizeof(bits));
    little_endian(out,bits,4);
  }

  void write_vtk(int t) const {
    std::vector<float> ux(n_,0),uy(n_,0),uz(n_,0),rho(n_,1),wmag(n_,0);
    const int nx=cfg_.nx,ny=cfg_.ny,nz=cfg_.nz;
    for (std::size_t i=0;i<n_;++i) if (!solid_[i]) {
      const auto m=macros(&f_[i*Q]);
      ux[i]=static_cast<float>(m.ux);
      uy[i]=static_cast<float>(m.uy);
      uz[i]=static_cast<float>(m.uz);
      rho[i]=static_cast<float>(m.rho);
    }
    // Three-dimensional curl: x/y periodic, z one-sided at inlet/outlet.
    for(int z=0;z<nz;++z) for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) {
      const auto i=idx(x,y,z,nx,ny);
      if (solid_[i]) continue;
      const auto xm=idx((x-1+nx)%nx,y,z,nx,ny);
      const auto xp=idx((x+1)%nx,y,z,nx,ny);
      const auto ym=idx(x,(y-1+ny)%ny,z,nx,ny);
      const auto yp=idx(x,(y+1)%ny,z,nx,ny);
      const auto zm=idx(x,y,std::max(z-1,0),nx,ny);
      const auto zp=idx(x,y,std::min(z+1,nz-1),nx,ny);
      const double dz=(z==0 || z==nz-1)?1.0:2.0;
      const double wx=(uz[yp]-uz[ym])/2-(uy[zp]-uy[zm])/dz;
      const double wy=(ux[zp]-ux[zm])/dz-(uz[xp]-uz[xm])/2;
      const double wz=(uy[xp]-uy[xm])/2-(ux[yp]-ux[ym])/2;
      wmag[i]=static_cast<float>(std::sqrt(wx*wx+wy*wy+wz*wz));
    }

    std::ostringstream name;
    name << "flow_" << std::setw(6) << std::setfill('0') << t << ".vti";
    const auto filename=(std::filesystem::path(cfg_.out)/name.str());
    std::ofstream out(filename, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write " + filename.string());

    // Offsets count bytes from the first byte after the appended-data underscore.
    // Each raw block starts with a UInt64 byte count, followed by its payload.
    const std::uint64_t count=static_cast<std::uint64_t>(n_);
    const std::uint64_t vector_bytes=count*3*4, scalar_bytes=count*4;
    const std::uint64_t speed_offset=8+vector_bytes;
    const std::uint64_t density_offset=speed_offset+8+scalar_bytes;
    const std::uint64_t curl_offset=density_offset+8+scalar_bytes;
    const std::uint64_t solid_offset=curl_offset+8+scalar_bytes;
    out << "<?xml version=\"1.0\"?>\n"
           "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
           "  <ImageData WholeExtent=\"0 " << nx-1 << " 0 " << ny-1 << " 0 " << nz-1
        << "\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n"
           "    <Piece Extent=\"0 " << nx-1 << " 0 " << ny-1 << " 0 " << nz-1 << "\">\n"
           "      <PointData Scalars=\"speed\" Vectors=\"velocity\">\n"
           "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"0\"/>\n"
           "        <DataArray type=\"Float32\" Name=\"speed\" format=\"appended\" offset=\"" << speed_offset << "\"/>\n"
           "        <DataArray type=\"Float32\" Name=\"density\" format=\"appended\" offset=\"" << density_offset << "\"/>\n"
           "        <DataArray type=\"Float32\" Name=\"vorticity_magnitude\" format=\"appended\" offset=\"" << curl_offset << "\"/>\n"
           "        <DataArray type=\"UInt8\" Name=\"solid\" format=\"appended\" offset=\"" << solid_offset << "\"/>\n"
           "      </PointData>\n      <CellData/>\n    </Piece>\n  </ImageData>\n"
           "  <AppendedData encoding=\"raw\">_";
    little_endian(out,vector_bytes,8);
    for(std::size_t i=0;i<n_;++i) {
      write_float(out,ux[i]); write_float(out,uy[i]); write_float(out,uz[i]);
    }
    little_endian(out,scalar_bytes,8);
    for(std::size_t i=0;i<n_;++i)
      write_float(out,std::sqrt(ux[i]*ux[i]+uy[i]*uy[i]+uz[i]*uz[i]));
    little_endian(out,scalar_bytes,8);
    for (float v:rho) write_float(out,v);
    little_endian(out,scalar_bytes,8);
    for (float v:wmag) write_float(out,v);
    little_endian(out,count,8);
    out.write(reinterpret_cast<const char*>(solid_.data()),
              static_cast<std::streamsize>(solid_.size()));
    out << "</AppendedData>\n</VTKFile>\n";
    if(!out) throw std::runtime_error("Error writing " + filename.string());
    std::cout << "wrote " << filename.string() << '\n';
  }
};
} // namespace

int main(int argc,char** argv) {
  try { Solver(parse_args(argc,argv)).run(); }
  catch(const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n'; return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
