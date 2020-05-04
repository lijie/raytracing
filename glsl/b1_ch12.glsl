// 啊... GLSL 写不了 interface
// 啊... 也没有 void* 或者 any 这样的类型
// 啊... 数组作为函数参数时必须定长
// 啊... 函数不能递归!!!
// :(

// 没有多态, 材质功能就只能靠type区分了
// 替换了一个看起来更好的随机函数
// 参考别人的代码学会了个小技巧, 可以把self设置为一张纹理, 估计是上一帧, 与当前帧算平均, 显示效果会更好一点

// 折射部分参考教程最新的代码总是不对
// 最后还是参考自己原来的C++代码写出来了...

#iChannel0 "self"

float pi = 3.14159265358;
float gSeed = 1.0;

// 射线
struct ray {
    vec3 origin;    // 原点
    vec3 direction; // 方向, 单位向量
};

ray ray_new(vec3 origin, vec3 direction) {
    ray r;
    r.origin = origin;
    r.direction = direction;
    return r;
}

vec3 ray_point_at(ray r, float t) {
    return r.origin + r.direction * t;
}

struct camera {
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 origin;
    vec3 u;
    vec3 v;
    vec3 w;
    float lens_radius;
};

camera camera_new(vec3 lookfrom, vec3 lookat, vec3 vup, float fov, float aspect, float aperture, float focus_dist) {
    camera c;

    c.origin = lookfrom;
    c.lens_radius = aperture / 2.0;

    float half_height = tan(radians(fov) / 2.0);
    float half_width = aspect * half_height;

    c.w = normalize(lookfrom - lookat);
    c.u = normalize(cross(vup, c.w));
    c.v = cross(c.w, c.u);

    c.lower_left_corner = c.origin - half_width * focus_dist * c.u - half_height * focus_dist * c.v - focus_dist * c.w;
    c.horizontal = 2.0 * half_width * focus_dist * c.u;
    c.vertical = 2.0 * half_height * focus_dist * c.v;
    return c;
}

vec2 randomInUnitDisk(inout float seed);

ray camera_get_ray(camera c, vec2 uv) {
    vec2 rd = c.lens_radius * randomInUnitDisk(gSeed);
    vec3 offset = c.u * rd.x + c.v * rd.y;
    vec3 direction = normalize(c.lower_left_corner + uv.x * c.horizontal + uv.y * c.vertical - c.origin - offset);
    return ray_new(c.origin + offset, direction);
}

int mat_type_lambertian = 1;
int mat_type_metal = 2;
int mat_type_dielectric = 3;
struct material {
    int type;
    float param1;
    vec3 albedo;
};

material material_new_lambertian(vec3 albedo) {
    material mat;
    mat.type = mat_type_lambertian;
    mat.albedo = albedo;
    return mat;
}

material material_new_metal(vec3 albedo, float fuzz) {
    material mat;
    mat.type = mat_type_metal;
    mat.albedo = albedo;
    mat.param1 = fuzz;
    return mat;
}

material material_new_dielectric(float ri) {
    material mat;
    mat.type = mat_type_dielectric;
    mat.param1 = ri;
    return mat;
}

struct hit_record {
    vec3 p;         // 交点
    vec3 normal;    // 交点处的法线
    float t;        // 交点与射线起点的距离
    bool front;    // 交点是在几何表面, 还是在里面
    material mat;
};

void hit_record_set_face_normal(out hit_record rec, ray r, vec3 out_normal) {
    // 如果在表面
    // 则射线(来自观察者)与法线的点积应该<0
    rec.front = dot(r.direction, out_normal) < 0.0;
    if (rec.front)
        rec.normal = out_normal;
    else
        rec.normal = -out_normal;
}

struct sphere {
    vec3 center;
    float radius;
    material mat;
};

sphere sphere_new(vec3 center, float radius, material mat) {
    sphere sp;
    sp.center = center;
    sp.radius = radius;
    sp.mat = mat;
    return sp;
}

bool sphere_hit(sphere sp, ray r, float t_min, float t_max, out hit_record rec) {
    // 计算一元二次方程
    // b^2 - 4ac
    float a = dot(r.direction, r.direction);
    float b = 2.0 * dot(r.origin - sp.center, r.direction);
    float c = dot((r.origin - sp.center), (r.origin - sp.center)) - sp.radius * sp.radius;
    float discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0)
        return false;
    
    float root = (-(b) - sqrt(discriminant)) / (2.0 * a);

    if (root <= t_min || root >= t_max) {
        root = (-(b) + sqrt(discriminant)) / (2.0 * a);
        if (root <= t_min || root >= t_max) {
            // 无解
            return false;
        }
    }

    // 有效解
    rec.t = root;
    rec.p = ray_point_at(r, rec.t);
    vec3 out_normal = (rec.p - sp.center) / sp.radius;
    // hit_record_set_face_normal(rec, r, out_normal);
    rec.normal = out_normal;
    rec.mat = sp.mat;
    return true;
}

sphere[100] world_sphere_array;
int world_sphere_count = 0;

bool world_hit(ray r, float t_min, float t_max, out hit_record rec) {
    float closest = t_max;
    bool hit = false;
    for (int i = 0; i < world_sphere_count; i++) {
        sphere sp = world_sphere_array[i];
        hit = sphere_hit(sp, r, t_min, closest, rec);
        if (hit) {
            closest = rec.t;
        }
    }
    return (t_max - closest) > 0.0;
}

void world_add(sphere sp) {
    world_sphere_array[world_sphere_count] = sp;
    world_sphere_count++;
}

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

bool material_scatter(material mat, ray r_in, hit_record rec, out vec3 attenuation, out ray scattered);

vec3 ray_color_depth(ray r, int depth) {
    hit_record rec;
    vec3 color = vec3(1.0);

    while (depth > 0) {
        depth--;

        if (!world_hit(r, 0.0001, 100000.0, rec))
            break;

        vec3 attenuation;
        ray scattered;
        if (material_scatter(rec.mat, r, rec, attenuation, scattered)) {
            color = color * attenuation;
            r = scattered;
        } else {
            return vec3(0);
        }
    }

    // scale t to [0, 1]
    float t = 0.5 * (r.direction.y + 1.0);
    return color * mix(vec3(1.0, 1.0, 1.0), vec3(0.5, 0.7, 1.0), t);
}

bool material_lambertian_scatter(material mat, ray r_in, hit_record rec, out vec3 attenuation, out ray scattered) {
    vec3 scatter_direction = normalize(rec.normal + randomInUnitSphere(gSeed));
    scattered = ray_new(rec.p, scatter_direction);
    attenuation = mat.albedo;
    return true;
}

// 反射公式
// vec3 reflect(vec3 v, vec3 n) {
//     return v - 2.0 * dot(v,n) * n;
// }

bool material_metal_scatter(material mat, ray r_in, out hit_record rec, out vec3 attenuation, out ray scattered) {
    vec3 reflected = reflect(r_in.direction, rec.normal);
    scattered = ray_new(rec.p, normalize(reflected + mat.param1 * randomInUnitSphere(gSeed)));
    attenuation = mat.albedo;
    return dot(scattered.direction, rec.normal) > 0.0;
}

float schlick(float cosine, float ref_idx) {
    float r0 = (1.0-ref_idx) / (1.0+ref_idx);
    r0 = r0*r0;
    return r0 + (1.0-r0)*pow((1.0 - cosine),5.0);
}


bool refract(vec3 v, vec3 n, float ni_over_nt, out vec3 refracted) {
  vec3 unit_v = normalize(v);
  float dt = dot(unit_v, n);
  float discriminat = 1.0 - (ni_over_nt * ni_over_nt * (1.0 - dt * dt));
  if (discriminat <= 0.0) return false;
  refracted = ni_over_nt * (unit_v - n * dt) / 1.0 - n * sqrt(discriminat);
  return true;
}

bool material_dielectric_scatter(material mat, ray r_in, out hit_record rec, out vec3 attenuation, out ray scattered) {
  vec3 outward_normal;
  vec3 reflected = reflect(r_in.direction, rec.normal);
  vec3 refracted;
  float ni_over_nt;
  float reflect_prob;
  float cosine;

  attenuation = vec3(1.0, 1.0, 1.0);  // 全反射

  // 这里我们总是假设材质外是空气， 折射率是1， 简化了计算
  if (dot(r_in.direction, rec.normal) > 0.0) {
    // 材质表面由外至内发生折射
    outward_normal = -rec.normal;
    ni_over_nt = mat.param1;
    // ni_over_nt = 1 / mat.param1;
    cosine = mat.param1 * dot(r_in.direction, rec.normal) / length(r_in.direction);
  } else {
    // 材质由内至外发生折射
    // 由于外部是空气，折射率是1， 根据 Snell's Law: 入射介质折射率*sin(入射角)
    // = 折射介质折射率*sin(折射角) 两边都除以入射介质折射率:
    outward_normal = rec.normal;
    // ni_over_nt = ref_idx_;
    ni_over_nt = 1.0 / mat.param1;
    cosine = -dot(r_in.direction, rec.normal) / length(r_in.direction);
  }

  // 计算折射射线
  if (refract(r_in.direction, outward_normal, ni_over_nt, refracted)) {
    reflect_prob = schlick(cosine, mat.param1);
    if (hash1(gSeed) >= reflect_prob) {
      scattered = ray_new(rec.p, normalize(refracted));
      return true;
    }
  }
  scattered = ray_new(rec.p, normalize(reflected));
  return true;
}

// vec3 refract(vec3 uv, vec3 n, float etai_over_etat) {
//     float cos_theta = dot(-uv, n);
//     vec3 r_out_parallel =  etai_over_etat * (uv + cos_theta*n);
//     vec3 r_out_perp = -sqrt(1.0 - pow(length(r_out_parallel),2.0)) * n;
//     return r_out_parallel + r_out_perp;
// }
// 
// bool material_dielectric_scatter(material mat, ray r_in, out hit_record rec, out vec3 attenuation, out ray scattered) {
//     attenuation = vec3(1.0, 1.0, 1.0);
//     float etai_over_etat = mat.param1;
//     if (rec.front) {
//         etai_over_etat = 1.0 / mat.param1;
//     }
//     vec3 unit_direction = r_in.direction;
//     float cos_theta = min(dot(-unit_direction, rec.normal), 1.0);
//     float sin_theta = sqrt(1.0 - cos_theta*cos_theta);
//     if (etai_over_etat * sin_theta > 1.0) {
//         vec3 reflected = reflect(unit_direction, rec.normal);
//         scattered = ray_new(rec.p, reflected);
//         return true;
//     }
// 
//     float reflect_prob = schlick(cos_theta, etai_over_etat);
//     if (hash1(gSeed) < 0.5)
//     {
//         vec3 reflected = reflect(unit_direction, rec.normal);
//         scattered = ray_new(rec.p, reflected);
//         return true;
//     }
// 
//     vec3 refracted = refract(unit_direction, rec.normal, etai_over_etat);
//     scattered = ray_new(rec.p, normalize(refracted));
//     return true;
// }

bool material_scatter(material mat, ray r_in, hit_record rec, out vec3 attenuation, out ray scattered) {
    if (mat.type == mat_type_lambertian) {
        return material_lambertian_scatter(mat, r_in, rec, attenuation, scattered);
    } else if (mat.type == mat_type_metal) {
        return material_metal_scatter(mat, r_in, rec, attenuation, scattered);
    } else if (mat.type == mat_type_dielectric) {
        return material_dielectric_scatter(mat, r_in, rec, attenuation, scattered);
    } else {
        return false;
    }
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

// 调高亮度
vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

void random_scene3() {
world_add(sphere_new(vec3(-10.354122, 0.200000, 3.016782), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(-10.595965, 0.200000, 7.650301), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(-8.545674, 0.200000, -8.576409), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(-8.192590, 0.200000, -7.492444), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(-8.204923, 0.200000, 4.790902), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(-7.727778, 0.200000, -2.989673), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(-6.318552, 0.200000, 3.376733), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(-5.135240, 0.200000, -7.218766), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(-5.423145, 0.200000, 7.877615), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(-3.982971, 0.200000, -10.516779), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(-2.729261, 0.200000, -9.811359), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(0.026148, 0.200000, 8.215833), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(2.755663, 0.200000, -9.457753), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(4.268294, 0.200000, -6.227857), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(4.879070, 0.200000, -2.747636), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(5.500525, 0.200000, -6.824818), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(5.095117, 0.200000, -5.611100), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(5.260274, 0.200000, -3.378430), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(5.710508, 0.200000, 9.770714), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(7.566857, 0.200000, 1.778295), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(7.337071, 0.200000, 9.452733), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(8.564330, 0.200000, -10.616593), 0.200000, material_new_dielectric(1.500000)));
// world_add(sphere_new(vec3(8.156615, 0.200000, -3.199265), 0.200000, material_new_dielectric(1.500000)));
world_add(sphere_new(vec3(9.479485, 0.200000, 6.441554), 0.200000, material_new_dielectric(1.500000)));
}

void random_scene2() {
world_add(sphere_new(vec3(-10.587423, 0.200000, -6.477859), 0.200000, material_new_metal(vec3(0.925657, 0.982940, 0.905164), 0.340556)));
world_add(sphere_new(vec3(-10.655239, 0.200000, -1.856652), 0.200000, material_new_metal(vec3(0.805048, 0.989776, 0.565401), 0.165136)));
// world_add(sphere_new(vec3(-10.859233, 0.200000, 1.819660), 0.200000, material_new_metal(vec3(0.803049, 0.542512, 0.938993), 0.172079)));
// world_add(sphere_new(vec3(-10.280071, 0.200000, 3.361077), 0.200000, material_new_metal(vec3(0.890988, 0.900754, 0.726066), 0.207160)));
// world_add(sphere_new(vec3(-10.581765, 0.200000, 7.083334), 0.200000, material_new_metal(vec3(0.501999, 0.540864, 0.648396), 0.272011)));
// world_add(sphere_new(vec3(-10.343547, 0.200000, 10.329765), 0.200000, material_new_metal(vec3(0.894848, 0.681295, 0.797281), 0.174795)));
// world_add(sphere_new(vec3(-9.837260, 0.200000, -6.667214), 0.200000, material_new_metal(vec3(0.525544, 0.880413, 0.648503), 0.139576)));
// world_add(sphere_new(vec3(-8.649608, 0.200000, -8.545042), 0.200000, material_new_metal(vec3(0.688482, 0.957732, 0.845790), 0.259102)));
// world_add(sphere_new(vec3(-8.114887, 0.200000, 9.455919), 0.200000, material_new_metal(vec3(0.511826, 0.751122, 0.503769), 0.411603)));
// world_add(sphere_new(vec3(-7.450420, 0.200000, -10.286883), 0.200000, material_new_metal(vec3(0.991562, 0.505921, 0.633747), 0.454573)));
// world_add(sphere_new(vec3(-7.813913, 0.200000, -8.267272), 0.200000, material_new_metal(vec3(0.722907, 0.880337, 0.526383), 0.412458)));
// world_add(sphere_new(vec3(-7.467559, 0.200000, 0.073006), 0.200000, material_new_metal(vec3(0.536256, 0.778832, 0.643239), 0.427381)));
// world_add(sphere_new(vec3(-7.973467, 0.200000, 4.470943), 0.200000, material_new_metal(vec3(0.940596, 0.666402, 0.530122), 0.128315)));
// world_add(sphere_new(vec3(-6.596130, 0.200000, -3.631178), 0.200000, material_new_metal(vec3(0.688787, 0.912458, 0.627445), 0.373730)));
// world_add(sphere_new(vec3(-5.580007, 0.200000, -4.929823), 0.200000, material_new_metal(vec3(0.971786, 0.539354, 0.556459), 0.067064)));
// world_add(sphere_new(vec3(-5.104395, 0.200000, 3.735282), 0.200000, material_new_metal(vec3(0.872997, 0.755379, 0.659062), 0.102969)));
world_add(sphere_new(vec3(-5.269167, 0.200000, 5.070479), 0.200000, material_new_metal(vec3(0.718101, 0.971496, 0.934309), 0.480651)));
world_add(sphere_new(vec3(-4.931388, 0.200000, -9.160234), 0.200000, material_new_metal(vec3(0.800317, 0.749474, 0.808344), 0.386242)));
// world_add(sphere_new(vec3(-4.489065, 0.200000, 4.391235), 0.200000, material_new_metal(vec3(0.865062, 0.753014, 0.957183), 0.439222)));
// world_add(sphere_new(vec3(-3.886068, 0.200000, -8.931416), 0.200000, material_new_metal(vec3(0.855312, 0.600055, 0.506836), 0.315775)));
// world_add(sphere_new(vec3(-3.187701, 0.200000, -3.408478), 0.200000, material_new_metal(vec3(0.501404, 0.822474, 0.752312), 0.118305)));
// world_add(sphere_new(vec3(-3.657518, 0.200000, -1.910239), 0.200000, material_new_metal(vec3(0.879101, 0.647526, 0.529847), 0.198874)));
// world_add(sphere_new(vec3(-2.980416, 0.200000, -4.959212), 0.200000, material_new_metal(vec3(0.668172, 0.950941, 0.865917), 0.081729)));
// world_add(sphere_new(vec3(-2.466268, 0.200000, -3.343190), 0.200000, material_new_metal(vec3(0.668157, 0.891552, 0.516739), 0.362606)));
// world_add(sphere_new(vec3(-2.373650, 0.200000, 3.772582), 0.200000, material_new_metal(vec3(0.745064, 0.960234, 0.522874), 0.400403)));
// world_add(sphere_new(vec3(-2.196930, 0.200000, 5.037959), 0.200000, material_new_metal(vec3(0.704932, 0.924207, 0.928938), 0.020905)));
// world_add(sphere_new(vec3(-1.496783, 0.200000, -4.656914), 0.200000, material_new_metal(vec3(0.711646, 0.878964, 0.508499), 0.283563)));
world_add(sphere_new(vec3(-1.251424, 0.200000, -3.703085), 0.200000, material_new_metal(vec3(0.917051, 0.725074, 0.610324), 0.090716)));
world_add(sphere_new(vec3(-0.414493, 0.200000, 3.189355), 0.200000, material_new_metal(vec3(0.634953, 0.733650, 0.532563), 0.247230)));
// world_add(sphere_new(vec3(-0.782519, 0.200000, 6.561831), 0.200000, material_new_metal(vec3(0.792962, 0.993912, 0.912961), 0.278558)));
// world_add(sphere_new(vec3(-0.284329, 0.200000, 8.044276), 0.200000, material_new_metal(vec3(0.686834, 0.537645, 0.690359), 0.312296)));
// world_add(sphere_new(vec3(0.727070, 0.200000, -7.762743), 0.200000, material_new_metal(vec3(0.786279, 0.525346, 0.990005), 0.409085)));
// world_add(sphere_new(vec3(0.495883, 0.200000, 8.704520), 0.200000, material_new_metal(vec3(0.777017, 0.708808, 0.576128), 0.213431)));
// world_add(sphere_new(vec3(1.663402, 0.200000, -8.650432), 0.200000, material_new_metal(vec3(0.515885, 0.927610, 0.669820), 0.376797)));
// world_add(sphere_new(vec3(1.857674, 0.200000, -5.693692), 0.200000, material_new_metal(vec3(0.708457, 0.625782, 0.701651), 0.183050)));
// world_add(sphere_new(vec3(1.137965, 0.200000, -1.817951), 0.200000, material_new_metal(vec3(0.968017, 0.718970, 0.875622), 0.134358)));
// world_add(sphere_new(vec3(1.699164, 0.200000, 2.525657), 0.200000, material_new_metal(vec3(0.665563, 0.985488, 0.810282), 0.174032)));
// world_add(sphere_new(vec3(2.665105, 0.200000, -9.529606), 0.200000, material_new_metal(vec3(0.802789, 0.627949, 0.859981), 0.000824)));
// world_add(sphere_new(vec3(2.259835, 0.200000, 5.345091), 0.200000, material_new_metal(vec3(0.641743, 0.938765, 0.737480), 0.115986)));
// world_add(sphere_new(vec3(2.528376, 0.200000, 9.774999), 0.200000, material_new_metal(vec3(0.757668, 0.819575, 0.558901), 0.361995)));
// world_add(sphere_new(vec3(4.418702, 0.200000, -8.288092), 0.200000, material_new_metal(vec3(0.984588, 0.843715, 0.871914), 0.467147)));
world_add(sphere_new(vec3(5.134559, 0.200000, -9.172100), 0.200000, material_new_metal(vec3(0.620838, 0.770211, 0.827006), 0.438063)));
world_add(sphere_new(vec3(5.860997, 0.200000, -6.871813), 0.200000, material_new_metal(vec3(0.684927, 0.643818, 0.824305), 0.096530)));
// world_add(sphere_new(vec3(5.386181, 0.200000, -5.721818), 0.200000, material_new_metal(vec3(0.801096, 0.958708, 0.858440), 0.377361)));
// world_add(sphere_new(vec3(5.862261, 0.200000, -2.930921), 0.200000, material_new_metal(vec3(0.628010, 0.966292, 0.857112), 0.266640)));
// world_add(sphere_new(vec3(5.538127, 0.200000, -1.274578), 0.200000, material_new_metal(vec3(0.947768, 0.714789, 0.803613), 0.373531)));
// world_add(sphere_new(vec3(5.890689, 0.200000, 3.123133), 0.200000, material_new_metal(vec3(0.938719, 0.712012, 0.970122), 0.284021)));
// world_add(sphere_new(vec3(6.181143, 0.200000, -5.126148), 0.200000, material_new_metal(vec3(0.616352, 0.765160, 0.756661), 0.263573)));
// world_add(sphere_new(vec3(6.293371, 0.200000, -3.719840), 0.200000, material_new_metal(vec3(0.975738, 0.813913, 0.609195), 0.194433)));
// world_add(sphere_new(vec3(6.032081, 0.200000, 6.350282), 0.200000, material_new_metal(vec3(0.583422, 0.848674, 0.507294), 0.287881)));
// world_add(sphere_new(vec3(6.640358, 0.200000, 10.494565), 0.200000, material_new_metal(vec3(0.977264, 0.989776, 0.876293), 0.221168)));
// world_add(sphere_new(vec3(7.821006, 0.200000, -5.361235), 0.200000, material_new_metal(vec3(0.630314, 0.782281, 0.774606), 0.441557)));
// world_add(sphere_new(vec3(7.544307, 0.200000, 5.089294), 0.200000, material_new_metal(vec3(0.612644, 0.780831, 0.512787), 0.458205)));
// world_add(sphere_new(vec3(7.206220, 0.200000, 6.660765), 0.200000, material_new_metal(vec3(0.747230, 0.770455, 0.597232), 0.497360)));
world_add(sphere_new(vec3(7.302463, 0.200000, 7.049358), 0.200000, material_new_metal(vec3(0.889950, 0.559267, 0.840373), 0.100574)));
world_add(sphere_new(vec3(8.618714, 0.200000, -6.272298), 0.200000, material_new_metal(vec3(0.977737, 0.880062, 0.702628), 0.367824)));
// world_add(sphere_new(vec3(10.585232, 0.200000, -5.609259), 0.200000, material_new_metal(vec3(0.588488, 0.500351, 0.945555), 0.314173)));
// world_add(sphere_new(vec3(10.411863, 0.200000, 3.334819), 0.200000, material_new_metal(vec3(0.884899, 0.997955, 0.689093), 0.385891)));
// world_add(sphere_new(vec3(10.611188, 0.200000, 9.050566), 0.200000, material_new_metal(vec3(0.508957, 0.758782, 0.829417), 0.188253)));
}

void random_scene1() {
world_add(sphere_new(vec3(-10.596405, 0.200000, -10.745000), 0.200000, material_new_lambertian(vec3(0.596459, 0.051637, 0.370106))));
world_add(sphere_new(vec3(-10.136613, 0.200000, -9.450667), 0.200000, material_new_lambertian(vec3(0.163987, 0.573336, 0.107150))));
world_add(sphere_new(vec3(-10.562511, 0.200000, -8.290701), 0.200000, material_new_lambertian(vec3(0.072692, 0.330116, 0.268589))));
world_add(sphere_new(vec3(-10.408231, 0.200000, -7.212064), 0.200000, material_new_lambertian(vec3(0.340794, 0.215936, 0.244234))));
world_add(sphere_new(vec3(-10.492444, 0.200000, -6.700888), 0.200000, material_new_lambertian(vec3(0.218700, 0.181551, 0.004897))));
world_add(sphere_new(vec3(-10.951906, 0.200000, -5.352776), 0.200000, material_new_lambertian(vec3(0.305900, 0.024606, 0.097403))));
world_add(sphere_new(vec3(-10.800510, 0.200000, -4.617307), 0.200000, material_new_lambertian(vec3(0.162408, 0.471710, 0.341639))));
world_add(sphere_new(vec3(-10.419987, 0.200000, -2.311850), 0.200000, material_new_lambertian(vec3(0.209423, 0.689729, 0.500788))));
world_add(sphere_new(vec3(-10.943364, 0.200000, -1.603958), 0.200000, material_new_lambertian(vec3(0.061045, 0.021502, 0.116194))));
world_add(sphere_new(vec3(-10.986486, 0.200000, -0.858602), 0.200000, material_new_lambertian(vec3(0.521261, 0.517877, 0.253979))));
world_add(sphere_new(vec3(-10.632221, 0.200000, 0.049303), 0.200000, material_new_lambertian(vec3(0.354702, 0.215027, 0.157666))));
world_add(sphere_new(vec3(-10.354506, 0.200000, 1.384039), 0.200000, material_new_lambertian(vec3(0.182230, 0.057946, 0.169206))));
// world_add(sphere_new(vec3(-10.262191, 0.200000, 2.023402), 0.200000, material_new_lambertian(vec3(0.477857, 0.042255, 0.227750))));
// world_add(sphere_new(vec3(-10.886975, 0.200000, 5.551970), 0.200000, material_new_lambertian(vec3(0.057510, 0.412187, 0.388834))));
// world_add(sphere_new(vec3(-10.666308, 0.200000, 7.670763), 0.200000, material_new_lambertian(vec3(0.007297, 0.034281, 0.479070))));
// world_add(sphere_new(vec3(-10.886370, 0.200000, 9.598773), 0.200000, material_new_lambertian(vec3(0.012903, 0.214934, 0.137591))));
// world_add(sphere_new(vec3(-9.861074, 0.200000, -10.663561), 0.200000, material_new_lambertian(vec3(0.051018, 0.201986, 0.575306))));
// world_add(sphere_new(vec3(-9.221320, 0.200000, -9.724921), 0.200000, material_new_lambertian(vec3(0.002786, 0.915170, 0.260970))));
// world_add(sphere_new(vec3(-9.283532, 0.200000, -8.856624), 0.200000, material_new_lambertian(vec3(0.572509, 0.523574, 0.077214))));
// world_add(sphere_new(vec3(-9.971874, 0.200000, -7.507578), 0.200000, material_new_lambertian(vec3(0.070925, 0.737850, 0.129764))));
// world_add(sphere_new(vec3(-9.706931, 0.200000, -6.271090), 0.200000, material_new_lambertian(vec3(0.074532, 0.681698, 0.074532))));
// world_add(sphere_new(vec3(-9.907630, 0.200000, -5.519059), 0.200000, material_new_lambertian(vec3(0.083434, 0.447839, 0.405904))));
// world_add(sphere_new(vec3(-9.617335, 0.200000, -3.245326), 0.200000, material_new_lambertian(vec3(0.344328, 0.008236, 0.110344))));
// world_add(sphere_new(vec3(-9.211789, 0.200000, -2.382055), 0.200000, material_new_lambertian(vec3(0.042061, 0.498302, 0.015692))));
// world_add(sphere_new(vec3(-9.704843, 0.200000, -1.107031), 0.200000, material_new_lambertian(vec3(0.104550, 0.332155, 0.013062))));
// world_add(sphere_new(vec3(-9.223270, 0.200000, -0.132163), 0.200000, material_new_lambertian(vec3(0.223774, 0.585762, 0.521953))));
// world_add(sphere_new(vec3(-9.260021, 0.200000, 1.628053), 0.200000, material_new_lambertian(vec3(0.441725, 0.162963, 0.018875))));
// world_add(sphere_new(vec3(-9.963744, 0.200000, 3.135136), 0.200000, material_new_lambertian(vec3(0.324958, 0.025213, 0.454133))));
// world_add(sphere_new(vec3(-9.156197, 0.200000, 6.702570), 0.200000, material_new_lambertian(vec3(0.548476, 0.631849, 0.762033))));
// world_add(sphere_new(vec3(-9.205664, 0.200000, 7.475063), 0.200000, material_new_lambertian(vec3(0.417245, 0.085318, 0.009062))));
// world_add(sphere_new(vec3(-9.412735, 0.200000, 8.749455), 0.200000, material_new_lambertian(vec3(0.388170, 0.440786, 0.571611))));
// world_add(sphere_new(vec3(-9.136201, 0.200000, 9.450481), 0.200000, material_new_lambertian(vec3(0.839744, 0.228004, 0.409476))));
// world_add(sphere_new(vec3(-9.864342, 0.200000, 10.878768), 0.200000, material_new_lambertian(vec3(0.026171, 0.046629, 0.090457))));
// world_add(sphere_new(vec3(-8.736183, 0.200000, -10.808118), 0.200000, material_new_lambertian(vec3(0.787271, 0.794996, 0.039263))));
// world_add(sphere_new(vec3(-8.289520, 0.200000, -9.625300), 0.200000, material_new_lambertian(vec3(0.117697, 0.192812, 0.143665))));
// world_add(sphere_new(vec3(-8.719813, 0.200000, -8.564901), 0.200000, material_new_lambertian(vec3(0.037060, 0.129806, 0.134446))));
// world_add(sphere_new(vec3(-8.417185, 0.200000, -7.677129), 0.200000, material_new_lambertian(vec3(0.124605, 0.043253, 0.106017))));
// world_add(sphere_new(vec3(-8.945808, 0.200000, -6.832481), 0.200000, material_new_lambertian(vec3(0.029047, 0.270203, 0.100249))));
// world_add(sphere_new(vec3(-8.932487, 0.200000, -5.366839), 0.200000, material_new_lambertian(vec3(0.000027, 0.095353, 0.151770))));
// world_add(sphere_new(vec3(-8.628321, 0.200000, -4.647465), 0.200000, material_new_lambertian(vec3(0.093424, 0.142405, 0.629956))));
// world_add(sphere_new(vec3(-8.193716, 0.200000, -2.953499), 0.200000, material_new_lambertian(vec3(0.259005, 0.084893, 0.773448))));
// world_add(sphere_new(vec3(-8.478298, 0.200000, 0.018375), 0.200000, material_new_lambertian(vec3(0.393175, 0.727391, 0.012558))));
// world_add(sphere_new(vec3(-8.650432, 0.200000, 1.543840), 0.200000, material_new_lambertian(vec3(0.316055, 0.213233, 0.126590))));
// world_add(sphere_new(vec3(-8.228626, 0.200000, 2.068337), 0.200000, material_new_lambertian(vec3(0.260523, 0.400846, 0.068247))));
// world_add(sphere_new(vec3(-8.675948, 0.200000, 3.266921), 0.200000, material_new_lambertian(vec3(0.785930, 0.000746, 0.249230))));
// world_add(sphere_new(vec3(-8.290866, 0.200000, 4.637254), 0.200000, material_new_lambertian(vec3(0.073725, 0.062281, 0.012241))));
// world_add(sphere_new(vec3(-8.774801, 0.200000, 5.732179), 0.200000, material_new_lambertian(vec3(0.060014, 0.380156, 0.492858))));
// world_add(sphere_new(vec3(-8.126533, 0.200000, 6.018842), 0.200000, material_new_lambertian(vec3(0.189132, 0.090630, 0.067211))));
// world_add(sphere_new(vec3(-8.877993, 0.200000, 8.135411), 0.200000, material_new_lambertian(vec3(0.295361, 0.025490, 0.088310))));
// world_add(sphere_new(vec3(-8.607282, 0.200000, 9.175457), 0.200000, material_new_lambertian(vec3(0.583163, 0.217841, 0.452663))));
// world_add(sphere_new(vec3(-8.351952, 0.200000, 10.184466), 0.200000, material_new_lambertian(vec3(0.343633, 0.157539, 0.127605))));
// world_add(sphere_new(vec3(-7.918589, 0.200000, -10.935591), 0.200000, material_new_lambertian(vec3(0.858295, 0.487796, 0.157161))));
// world_add(sphere_new(vec3(-7.133125, 0.200000, -9.329099), 0.200000, material_new_lambertian(vec3(0.272948, 0.232673, 0.506314))));
// world_add(sphere_new(vec3(-7.587918, 0.200000, -7.113926), 0.200000, material_new_lambertian(vec3(0.153257, 0.572908, 0.176254))));
// world_add(sphere_new(vec3(-7.311576, 0.200000, -6.908701), 0.200000, material_new_lambertian(vec3(0.349751, 0.168144, 0.140295))));
// world_add(sphere_new(vec3(-7.826548, 0.200000, -5.809024), 0.200000, material_new_lambertian(vec3(0.208450, 0.231795, 0.160976))));
// world_add(sphere_new(vec3(-7.802542, 0.200000, -4.625300), 0.200000, material_new_lambertian(vec3(0.018413, 0.224216, 0.126012))));
// world_add(sphere_new(vec3(-7.510654, 0.200000, -3.761754), 0.200000, material_new_lambertian(vec3(0.502012, 0.201243, 0.585322))));
// world_add(sphere_new(vec3(-7.281884, 0.200000, -0.435038), 0.200000, material_new_lambertian(vec3(0.558817, 0.097500, 0.135488))));
// world_add(sphere_new(vec3(-7.870907, 0.200000, 1.681420), 0.200000, material_new_lambertian(vec3(0.049590, 0.047482, 0.253390))));
// world_add(sphere_new(vec3(-7.372332, 0.200000, 4.481188), 0.200000, material_new_lambertian(vec3(0.045969, 0.131603, 0.086966))));
// world_add(sphere_new(vec3(-7.507111, 0.200000, 6.246294), 0.200000, material_new_lambertian(vec3(0.746619, 0.390407, 0.081181))));
// world_add(sphere_new(vec3(-7.568471, 0.200000, 7.856246), 0.200000, material_new_lambertian(vec3(0.524943, 0.351239, 0.038248))));
// world_add(sphere_new(vec3(-7.316272, 0.200000, 9.113767), 0.200000, material_new_lambertian(vec3(0.269157, 0.187184, 0.428471))));
// world_add(sphere_new(vec3(-6.568966, 0.200000, -10.447151), 0.200000, material_new_lambertian(vec3(0.290641, 0.147754, 0.225751))));
// world_add(sphere_new(vec3(-6.867226, 0.200000, -9.625822), 0.200000, material_new_lambertian(vec3(0.100718, 0.179011, 0.567477))));
// world_add(sphere_new(vec3(-6.940123, 0.200000, -8.155318), 0.200000, material_new_lambertian(vec3(0.047680, 0.490056, 0.017341))));
// world_add(sphere_new(vec3(-6.425562, 0.200000, -7.741704), 0.200000, material_new_lambertian(vec3(0.359032, 0.152693, 0.127426))));
// world_add(sphere_new(vec3(-6.310038, 0.200000, -6.170122), 0.200000, material_new_lambertian(vec3(0.484287, 0.016151, 0.013410))));
// world_add(sphere_new(vec3(-6.704816, 0.200000, -5.953087), 0.200000, material_new_lambertian(vec3(0.230846, 0.035611, 0.221531))));
// world_add(sphere_new(vec3(-6.442509, 0.200000, -4.886728), 0.200000, material_new_lambertian(vec3(0.276219, 0.449750, 0.406981))));
// world_add(sphere_new(vec3(-6.999066, 0.200000, -3.606815), 0.200000, material_new_lambertian(vec3(0.258688, 0.292160, 0.134634))));
// world_add(sphere_new(vec3(-6.507111, 0.200000, -0.335499), 0.200000, material_new_lambertian(vec3(0.216314, 0.192819, 0.903510))));
// world_add(sphere_new(vec3(-6.583248, 0.200000, 0.437873), 0.200000, material_new_lambertian(vec3(0.050980, 0.176993, 0.240632))));
// world_add(sphere_new(vec3(-6.180258, 0.200000, 1.196112), 0.200000, material_new_lambertian(vec3(0.543886, 0.376768, 0.843023))));
// world_add(sphere_new(vec3(-6.477557, 0.200000, 2.086657), 0.200000, material_new_lambertian(vec3(0.148109, 0.392794, 0.301686))));
// world_add(sphere_new(vec3(-6.619944, 0.200000, 4.738936), 0.200000, material_new_lambertian(vec3(0.226541, 0.401206, 0.562521))));
// world_add(sphere_new(vec3(-6.477172, 0.200000, 5.571966), 0.200000, material_new_lambertian(vec3(0.436999, 0.288397, 0.098532))));
// world_add(sphere_new(vec3(-6.803998, 0.200000, 6.153594), 0.200000, material_new_lambertian(vec3(0.157424, 0.155062, 0.366437))));
// world_add(sphere_new(vec3(-6.437236, 0.200000, 7.620637), 0.200000, material_new_lambertian(vec3(0.157050, 0.717944, 0.571566))));
// world_add(sphere_new(vec3(-6.408478, 0.200000, 8.744841), 0.200000, material_new_lambertian(vec3(0.274066, 0.100721, 0.425530))));
// world_add(sphere_new(vec3(-6.728849, 0.200000, 9.868770), 0.200000, material_new_lambertian(vec3(0.070274, 0.125309, 0.328662))));
// world_add(sphere_new(vec3(-5.165371, 0.200000, -10.197397), 0.200000, material_new_lambertian(vec3(0.089907, 0.737436, 0.118372))));
// world_add(sphere_new(vec3(-5.837645, 0.200000, -9.476733), 0.200000, material_new_lambertian(vec3(0.068665, 0.120370, 0.442533))));
// world_add(sphere_new(vec3(-5.929493, 0.200000, -8.322810), 0.200000, material_new_lambertian(vec3(0.077761, 0.510890, 0.493573))));
// world_add(sphere_new(vec3(-5.826960, 0.200000, -7.469619), 0.200000, material_new_lambertian(vec3(0.659315, 0.833108, 0.009993))));
// world_add(sphere_new(vec3(-5.661968, 0.200000, -6.891726), 0.200000, material_new_lambertian(vec3(0.041569, 0.093097, 0.078928))));
// world_add(sphere_new(vec3(-5.916199, 0.200000, -5.503787), 0.200000, material_new_lambertian(vec3(0.163449, 0.031603, 0.003955))));
world_add(sphere_new(vec3(-5.361318, 0.200000, -4.512906), 0.200000, material_new_lambertian(vec3(0.295903, 0.643802, 0.092613))));
world_add(sphere_new(vec3(-5.622883, 0.200000, -2.326463), 0.200000, material_new_lambertian(vec3(0.162860, 0.683393, 0.279259))));
world_add(sphere_new(vec3(-5.167403, 0.200000, -1.718302), 0.200000, material_new_lambertian(vec3(0.265093, 0.047107, 0.209495))));
world_add(sphere_new(vec3(-5.349425, 0.200000, -0.990167), 0.200000, material_new_lambertian(vec3(0.050025, 0.006671, 0.098710))));
world_add(sphere_new(vec3(-5.277709, 0.200000, 0.894369), 0.200000, material_new_lambertian(vec3(0.160189, 0.074722, 0.364680))));
world_add(sphere_new(vec3(-5.761480, 0.200000, 1.851274), 0.200000, material_new_lambertian(vec3(0.105479, 0.192615, 0.023910))));
world_add(sphere_new(vec3(-5.697674, 0.200000, 2.054439), 0.200000, material_new_lambertian(vec3(0.768354, 0.324915, 0.079213))));
world_add(sphere_new(vec3(-5.984069, 0.200000, 3.769341), 0.200000, material_new_lambertian(vec3(0.061653, 0.245995, 0.031115))));
// world_add(sphere_new(vec3(-5.210910, 0.200000, 4.367833), 0.200000, material_new_lambertian(vec3(0.422019, 0.224499, 0.029046))));
// world_add(sphere_new(vec3(-5.652190, 0.200000, 5.573614), 0.200000, material_new_lambertian(vec3(0.342112, 0.304984, 0.092067))));
// world_add(sphere_new(vec3(-5.662709, 0.200000, 6.447569), 0.200000, material_new_lambertian(vec3(0.075745, 0.058589, 0.509764))));
// world_add(sphere_new(vec3(-5.467696, 0.200000, 8.467125), 0.200000, material_new_lambertian(vec3(0.252466, 0.161005, 0.124070))));
// world_add(sphere_new(vec3(-5.703305, 0.200000, 9.130934), 0.200000, material_new_lambertian(vec3(0.365379, 0.967273, 0.768584))));
// world_add(sphere_new(vec3(-5.616044, 0.200000, 10.442955), 0.200000, material_new_lambertian(vec3(0.038083, 0.018587, 0.005073))));
// world_add(sphere_new(vec3(-4.439378, 0.200000, -9.990030), 0.200000, material_new_lambertian(vec3(0.599928, 0.279701, 0.261769))));
// world_add(sphere_new(vec3(-4.500766, 0.200000, -8.550865), 0.200000, material_new_lambertian(vec3(0.095857, 0.020861, 0.245697))));
// world_add(sphere_new(vec3(-4.919083, 0.200000, -7.248018), 0.200000, material_new_lambertian(vec3(0.118863, 0.269205, 0.503833))));
// world_add(sphere_new(vec3(-4.667818, 0.200000, -6.706189), 0.200000, material_new_lambertian(vec3(0.219354, 0.576057, 0.688203))));
// world_add(sphere_new(vec3(-4.487472, 0.200000, -5.538862), 0.200000, material_new_lambertian(vec3(0.251812, 0.670066, 0.226407))));
// world_add(sphere_new(vec3(-4.379803, 0.200000, -4.426167), 0.200000, material_new_lambertian(vec3(0.723401, 0.668188, 0.464274))));
// world_add(sphere_new(vec3(-4.850829, 0.200000, -3.796225), 0.200000, material_new_lambertian(vec3(0.154306, 0.327151, 0.044924))));
// world_add(sphere_new(vec3(-4.758321, 0.200000, -1.613242), 0.200000, material_new_lambertian(vec3(0.273681, 0.092605, 0.314200))));
// world_add(sphere_new(vec3(-4.819709, 0.200000, -0.919550), 0.200000, material_new_lambertian(vec3(0.040238, 0.198796, 0.228710))));
// world_add(sphere_new(vec3(-4.196243, 0.200000, 0.392938), 0.200000, material_new_lambertian(vec3(0.049289, 0.065448, 0.847507))));
// world_add(sphere_new(vec3(-4.497497, 0.200000, 1.866738), 0.200000, material_new_lambertian(vec3(0.091138, 0.736915, 0.576757))));
// world_add(sphere_new(vec3(-4.928834, 0.200000, 3.476217), 0.200000, material_new_lambertian(vec3(0.093523, 0.307867, 0.121106))));
// world_add(sphere_new(vec3(-4.821741, 0.200000, 4.584133), 0.200000, material_new_lambertian(vec3(0.725909, 0.087872, 0.000774))));
// world_add(sphere_new(vec3(-4.584292, 0.200000, 6.770879), 0.200000, material_new_lambertian(vec3(0.364229, 0.365820, 0.142879))));
// world_add(sphere_new(vec3(-4.358818, 0.200000, 7.237559), 0.200000, material_new_lambertian(vec3(0.468950, 0.002506, 0.214723))));
// world_add(sphere_new(vec3(-4.556413, 0.200000, 8.282522), 0.200000, material_new_lambertian(vec3(0.007625, 0.441832, 0.102027))));
// world_add(sphere_new(vec3(-4.276308, 0.200000, 9.098083), 0.200000, material_new_lambertian(vec3(0.431068, 0.343333, 0.211382))));
// world_add(sphere_new(vec3(-4.940754, 0.200000, 10.747835), 0.200000, material_new_lambertian(vec3(0.163762, 0.293749, 0.411629))));
// world_add(sphere_new(vec3(-3.729563, 0.200000, -9.855443), 0.200000, material_new_lambertian(vec3(0.752762, 0.582147, 0.597993))));
// world_add(sphere_new(vec3(-3.101593, 0.200000, -8.589016), 0.200000, material_new_lambertian(vec3(0.585067, 0.173299, 0.167727))));
// world_add(sphere_new(vec3(-3.940672, 0.200000, -7.795703), 0.200000, material_new_lambertian(vec3(0.485441, 0.494789, 0.376867))));
// world_add(sphere_new(vec3(-3.909525, 0.200000, -6.139689), 0.200000, material_new_lambertian(vec3(0.054784, 0.957325, 0.306872))));
// world_add(sphere_new(vec3(-3.246947, 0.200000, -5.997446), 0.200000, material_new_lambertian(vec3(0.668173, 0.437433, 0.031089))));
// world_add(sphere_new(vec3(-3.168007, 0.200000, -4.951878), 0.200000, material_new_lambertian(vec3(0.169823, 0.157955, 0.242096))));
// world_add(sphere_new(vec3(-3.928147, 0.200000, -3.702399), 0.200000, material_new_lambertian(vec3(0.009077, 0.297212, 0.076232))));
// world_add(sphere_new(vec3(-3.119474, 0.200000, -2.825312), 0.200000, material_new_lambertian(vec3(0.049639, 0.092068, 0.703438))));
// world_add(sphere_new(vec3(-3.527271, 0.200000, -1.100632), 0.200000, material_new_lambertian(vec3(0.169597, 0.135385, 0.357919))));
// world_add(sphere_new(vec3(-3.674657, 0.200000, -0.687896), 0.200000, material_new_lambertian(vec3(0.140933, 0.323284, 0.030421))));
// world_add(sphere_new(vec3(-3.905762, 0.200000, 1.881542), 0.200000, material_new_lambertian(vec3(0.024224, 0.476856, 0.118678))));
// world_add(sphere_new(vec3(-3.145100, 0.200000, 3.099210), 0.200000, material_new_lambertian(vec3(0.557505, 0.013831, 0.314314))));
// world_add(sphere_new(vec3(-3.665044, 0.200000, 4.563561), 0.200000, material_new_lambertian(vec3(0.731338, 0.262160, 0.695088))));
// world_add(sphere_new(vec3(-3.435807, 0.200000, 5.638325), 0.200000, material_new_lambertian(vec3(0.116404, 0.171830, 0.821248))));
// world_add(sphere_new(vec3(-3.244282, 0.200000, 6.495471), 0.200000, material_new_lambertian(vec3(0.117601, 0.867151, 0.475692))));
// world_add(sphere_new(vec3(-3.138426, 0.200000, 7.470312), 0.200000, material_new_lambertian(vec3(0.665911, 0.066176, 0.045062))));
// world_add(sphere_new(vec3(-3.791473, 0.200000, 8.373739), 0.200000, material_new_lambertian(vec3(0.094775, 0.201671, 0.617411))));
// world_add(sphere_new(vec3(-3.394964, 0.200000, 9.060894), 0.200000, material_new_lambertian(vec3(0.095202, 0.324741, 0.244924))));
// world_add(sphere_new(vec3(-3.625684, 0.200000, 10.446251), 0.200000, material_new_lambertian(vec3(0.407500, 0.091287, 0.870452))));
// world_add(sphere_new(vec3(-2.986486, 0.200000, -10.759777), 0.200000, material_new_lambertian(vec3(0.299050, 0.029194, 0.025501))));
// world_add(sphere_new(vec3(-2.195529, 0.200000, -9.292596), 0.200000, material_new_lambertian(vec3(0.055524, 0.310043, 0.007181))));
// world_add(sphere_new(vec3(-2.275347, 0.200000, -8.129005), 0.200000, material_new_lambertian(vec3(0.242854, 0.368582, 0.019207))));
// world_add(sphere_new(vec3(-2.355770, 0.200000, -4.445503), 0.200000, material_new_lambertian(vec3(0.277548, 0.068247, 0.025734))));
// world_add(sphere_new(vec3(-2.946110, 0.200000, -3.706848), 0.200000, material_new_lambertian(vec3(0.709661, 0.303218, 0.006460))));
// world_add(sphere_new(vec3(-2.819242, 0.200000, -2.973220), 0.200000, material_new_lambertian(vec3(0.575981, 0.257151, 0.384387))));
// world_add(sphere_new(vec3(-2.426441, 0.200000, -1.984701), 0.200000, material_new_lambertian(vec3(0.291836, 0.093880, 0.055446))));
// world_add(sphere_new(vec3(-2.214536, 0.200000, -0.131065), 0.200000, material_new_lambertian(vec3(0.098970, 0.087562, 0.414191))));
// world_add(sphere_new(vec3(-2.634474, 0.200000, 0.162932), 0.200000, material_new_lambertian(vec3(0.216136, 0.184911, 0.425348))));
// world_add(sphere_new(vec3(-2.285757, 0.200000, 1.065893), 0.200000, material_new_lambertian(vec3(0.175794, 0.540925, 0.135936))));
// world_add(sphere_new(vec3(-2.477639, 0.200000, 2.315729), 0.200000, material_new_lambertian(vec3(0.044513, 0.131719, 0.607052))));
// world_add(sphere_new(vec3(-2.980004, 0.200000, 3.043617), 0.200000, material_new_lambertian(vec3(0.346255, 0.195340, 0.003774))));
// world_add(sphere_new(vec3(-2.416965, 0.200000, 4.404859), 0.200000, material_new_lambertian(vec3(0.598560, 0.000046, 0.273349))));
// world_add(sphere_new(vec3(-2.557595, 0.200000, 5.343031), 0.200000, material_new_lambertian(vec3(0.074497, 0.205334, 0.291197))));
// world_add(sphere_new(vec3(-2.220661, 0.200000, 6.301639), 0.200000, material_new_lambertian(vec3(0.123486, 0.010454, 0.122968))));
// world_add(sphere_new(vec3(-2.573360, 0.200000, 7.518104), 0.200000, material_new_lambertian(vec3(0.156360, 0.086331, 0.440055))));
// world_add(sphere_new(vec3(-2.899939, 0.200000, 8.800323), 0.200000, material_new_lambertian(vec3(0.314222, 0.446732, 0.042864))));
// world_add(sphere_new(vec3(-2.287817, 0.200000, 9.110334), 0.200000, material_new_lambertian(vec3(0.213379, 0.313751, 0.371116))));
// world_add(sphere_new(vec3(-2.560808, 0.200000, 10.172848), 0.200000, material_new_lambertian(vec3(0.002978, 0.627435, 0.172630))));
// world_add(sphere_new(vec3(-1.618241, 0.200000, -10.996567), 0.200000, material_new_lambertian(vec3(0.252793, 0.297602, 0.057306))));
// world_add(sphere_new(vec3(-1.236949, 0.200000, -7.210718), 0.200000, material_new_lambertian(vec3(0.635166, 0.661278, 0.512847))));
// world_add(sphere_new(vec3(-1.310889, 0.200000, -6.663012), 0.200000, material_new_lambertian(vec3(0.377874, 0.003596, 0.062600))));
// world_add(sphere_new(vec3(-1.332560, 0.200000, -5.682870), 0.200000, material_new_lambertian(vec3(0.048601, 0.110887, 0.015502))));
// world_add(sphere_new(vec3(-1.969677, 0.200000, -4.388757), 0.200000, material_new_lambertian(vec3(0.556945, 0.183882, 0.203810))));
// world_add(sphere_new(vec3(-1.127934, 0.200000, -3.866594), 0.200000, material_new_lambertian(vec3(0.004912, 0.554241, 0.356217))));
// world_add(sphere_new(vec3(-1.927131, 0.200000, -2.436357), 0.200000, material_new_lambertian(vec3(0.006744, 0.407050, 0.501305))));
// world_add(sphere_new(vec3(-1.129279, 0.200000, -1.800455), 0.200000, material_new_lambertian(vec3(0.688937, 0.615795, 0.630964))));
// world_add(sphere_new(vec3(-1.155455, 0.200000, -0.819681), 0.200000, material_new_lambertian(vec3(0.070538, 0.078463, 0.746922))));
// world_add(sphere_new(vec3(-1.224781, 0.200000, 1.728883), 0.200000, material_new_lambertian(vec3(0.031330, 0.845398, 0.126026))));
// world_add(sphere_new(vec3(-1.202341, 0.200000, 2.696637), 0.200000, material_new_lambertian(vec3(0.027919, 0.654789, 0.109506))));
// world_add(sphere_new(vec3(-1.351567, 0.200000, 3.410572), 0.200000, material_new_lambertian(vec3(0.120997, 0.287127, 0.315357))));
// world_add(sphere_new(vec3(-1.985525, 0.200000, 4.809992), 0.200000, material_new_lambertian(vec3(0.125808, 0.098521, 0.314126))));
// world_add(sphere_new(vec3(-1.376122, 0.200000, 5.745994), 0.200000, material_new_lambertian(vec3(0.209210, 0.269547, 0.453526))));
// world_add(sphere_new(vec3(-1.782711, 0.200000, 6.895248), 0.200000, material_new_lambertian(vec3(0.230396, 0.496079, 0.306650))));
// world_add(sphere_new(vec3(-1.495712, 0.200000, 8.410489), 0.200000, material_new_lambertian(vec3(0.205576, 0.581928, 0.010093))));
// world_add(sphere_new(vec3(-1.743791, 0.200000, 10.825813), 0.200000, material_new_lambertian(vec3(0.426900, 0.325275, 0.580909))));
// world_add(sphere_new(vec3(-0.532325, 0.200000, -10.532051), 0.200000, material_new_lambertian(vec3(0.224648, 0.215735, 0.076780))));
// world_add(sphere_new(vec3(-0.963497, 0.200000, -7.186987), 0.200000, material_new_lambertian(vec3(0.108793, 0.481855, 0.288949))));
// world_add(sphere_new(vec3(-0.846791, 0.200000, -6.171413), 0.200000, material_new_lambertian(vec3(0.098922, 0.343607, 0.322241))));
// world_add(sphere_new(vec3(-0.962838, 0.200000, -3.799521), 0.200000, material_new_lambertian(vec3(0.211294, 0.301193, 0.548137))));
// world_add(sphere_new(vec3(-0.402326, 0.200000, -2.649937), 0.200000, material_new_lambertian(vec3(0.984595, 0.371507, 0.095166))));
// world_add(sphere_new(vec3(-0.727970, 0.200000, -0.675729), 0.200000, material_new_lambertian(vec3(0.300973, 0.002022, 0.046462))));
// world_add(sphere_new(vec3(-0.462917, 0.200000, 0.676366), 0.200000, material_new_lambertian(vec3(0.126234, 0.579318, 0.268532))));
// world_add(sphere_new(vec3(-0.592724, 0.200000, 1.087811), 0.200000, material_new_lambertian(vec3(0.072471, 0.112868, 0.129988))));
// world_add(sphere_new(vec3(-0.428968, 0.200000, 2.115827), 0.200000, material_new_lambertian(vec3(0.093369, 0.902707, 0.061495))));
// world_add(sphere_new(vec3(-0.293036, 0.200000, 3.283566), 0.200000, material_new_lambertian(vec3(0.273793, 0.465898, 0.106803))));
// world_add(sphere_new(vec3(-0.569625, 0.200000, 4.836222), 0.200000, material_new_lambertian(vec3(0.187985, 0.001335, 0.147547))));
// world_add(sphere_new(vec3(-0.477502, 0.200000, 5.012168), 0.200000, material_new_lambertian(vec3(0.251922, 0.362058, 0.232222))));
// world_add(sphere_new(vec3(-0.395926, 0.200000, 7.450591), 0.200000, material_new_lambertian(vec3(0.241456, 0.439408, 0.237904))));
// world_add(sphere_new(vec3(-0.707974, 0.200000, 8.633244), 0.200000, material_new_lambertian(vec3(0.605980, 0.011381, 0.253383))));
// world_add(sphere_new(vec3(-0.187509, 0.200000, 10.057900), 0.200000, material_new_lambertian(vec3(0.413805, 0.229793, 0.323581))));
world_add(sphere_new(vec3(0.308918, 0.200000, -10.500986), 0.200000, material_new_lambertian(vec3(0.016361, 0.297358, 0.719631))));
world_add(sphere_new(vec3(0.266646, 0.200000, -8.427073), 0.200000, material_new_lambertian(vec3(0.542022, 0.154104, 0.020685))));
world_add(sphere_new(vec3(0.382748, 0.200000, -7.869973), 0.200000, material_new_lambertian(vec3(0.176437, 0.190915, 0.149526))));
world_add(sphere_new(vec3(0.507886, 0.200000, -6.222117), 0.200000, material_new_lambertian(vec3(0.381329, 0.272556, 0.250092))));
world_add(sphere_new(vec3(0.287109, 0.200000, -4.947401), 0.200000, material_new_lambertian(vec3(0.398376, 0.034281, 0.518466))));
world_add(sphere_new(vec3(0.084625, 0.200000, -2.544438), 0.200000, material_new_lambertian(vec3(0.066106, 0.496242, 0.596243))));
world_add(sphere_new(vec3(0.153621, 0.200000, -1.697977), 0.200000, material_new_lambertian(vec3(0.168855, 0.024731, 0.591176))));
world_add(sphere_new(vec3(0.049303, 0.200000, -0.860195), 0.200000, material_new_lambertian(vec3(0.256016, 0.265978, 0.014602))));
world_add(sphere_new(vec3(0.673867, 0.200000, 0.650850), 0.200000, material_new_lambertian(vec3(0.052360, 0.000230, 0.328555))));
world_add(sphere_new(vec3(0.653597, 0.200000, 2.164882), 0.200000, material_new_lambertian(vec3(0.293640, 0.025781, 0.481442))));
// world_add(sphere_new(vec3(0.853801, 0.200000, 4.075396), 0.200000, material_new_lambertian(vec3(0.005569, 0.620463, 0.386989))));
// world_add(sphere_new(vec3(0.764946, 0.200000, 5.513901), 0.200000, material_new_lambertian(vec3(0.439547, 0.088323, 0.289049))));
// world_add(sphere_new(vec3(0.090695, 0.200000, 6.017634), 0.200000, material_new_lambertian(vec3(0.574473, 0.228976, 0.153610))));
// world_add(sphere_new(vec3(0.397525, 0.200000, 8.330891), 0.200000, material_new_lambertian(vec3(0.260698, 0.037154, 0.090470))));
// world_add(sphere_new(vec3(0.457018, 0.200000, 9.277111), 0.200000, material_new_lambertian(vec3(0.209209, 0.452274, 0.373237))));
// world_add(sphere_new(vec3(0.645109, 0.200000, 10.527140), 0.200000, material_new_lambertian(vec3(0.534403, 0.499016, 0.066597))));
// world_add(sphere_new(vec3(1.457210, 0.200000, -9.783892), 0.200000, material_new_lambertian(vec3(0.502080, 0.082150, 0.004610))));
// world_add(sphere_new(vec3(1.831196, 0.200000, -8.178115), 0.200000, material_new_lambertian(vec3(0.065509, 0.437162, 0.094666))));
// world_add(sphere_new(vec3(1.835920, 0.200000, -7.654552), 0.200000, material_new_lambertian(vec3(0.181007, 0.141105, 0.615169))));
// world_add(sphere_new(vec3(1.183972, 0.200000, -5.728135), 0.200000, material_new_lambertian(vec3(0.668247, 0.143537, 0.477920))));
// world_add(sphere_new(vec3(1.159059, 0.200000, -4.976791), 0.200000, material_new_lambertian(vec3(0.163404, 0.369309, 0.160290))));
// world_add(sphere_new(vec3(1.675845, 0.200000, -3.910733), 0.200000, material_new_lambertian(vec3(0.907158, 0.412654, 0.566955))));
// world_add(sphere_new(vec3(1.719434, 0.200000, -2.772549), 0.200000, material_new_lambertian(vec3(0.041384, 0.578658, 0.011237))));
// world_add(sphere_new(vec3(1.077017, 0.200000, -1.229423), 0.200000, material_new_lambertian(vec3(0.001354, 0.669453, 0.096962))));
// world_add(sphere_new(vec3(1.455699, 0.200000, 0.022605), 0.200000, material_new_lambertian(vec3(0.293328, 0.075036, 0.716288))));
// world_add(sphere_new(vec3(1.208664, 0.200000, 1.331935), 0.200000, material_new_lambertian(vec3(0.563070, 0.282820, 0.482114))));
// world_add(sphere_new(vec3(1.585974, 0.200000, 3.275573), 0.200000, material_new_lambertian(vec3(0.269872, 0.168857, 0.551908))));
// world_add(sphere_new(vec3(1.883630, 0.200000, 4.059026), 0.200000, material_new_lambertian(vec3(0.334208, 0.034920, 0.023870))));
// world_add(sphere_new(vec3(1.757558, 0.200000, 5.655107), 0.200000, material_new_lambertian(vec3(0.138973, 0.080437, 0.214846))));
// world_add(sphere_new(vec3(1.448668, 0.200000, 7.195315), 0.200000, material_new_lambertian(vec3(0.138080, 0.028593, 0.140659))));
// world_add(sphere_new(vec3(1.216355, 0.200000, 8.654805), 0.200000, material_new_lambertian(vec3(0.188153, 0.488207, 0.002355))));
// world_add(sphere_new(vec3(1.055812, 0.200000, 10.385961), 0.200000, material_new_lambertian(vec3(0.301384, 0.590077, 0.017809))));
// world_add(sphere_new(vec3(2.485446, 0.200000, -10.750850), 0.200000, material_new_lambertian(vec3(0.338271, 0.451711, 0.024562))));
// world_add(sphere_new(vec3(2.235472, 0.200000, -9.844676), 0.200000, material_new_lambertian(vec3(0.573427, 0.203801, 0.012424))));
// world_add(sphere_new(vec3(2.761403, 0.200000, -8.213904), 0.200000, material_new_lambertian(vec3(0.439673, 0.587016, 0.521064))));
// world_add(sphere_new(vec3(2.492642, 0.200000, -7.485824), 0.200000, material_new_lambertian(vec3(0.056776, 0.132421, 0.031982))));
// world_add(sphere_new(vec3(2.077676, 0.200000, -6.448222), 0.200000, material_new_lambertian(vec3(0.324073, 0.064984, 0.508408))));
// world_add(sphere_new(vec3(2.723911, 0.200000, -5.290371), 0.200000, material_new_lambertian(vec3(0.021000, 0.003839, 0.122416))));
// world_add(sphere_new(vec3(2.761678, 0.200000, -4.367360), 0.200000, material_new_lambertian(vec3(0.660145, 0.581064, 0.798190))));
// world_add(sphere_new(vec3(2.408786, 0.200000, -3.138591), 0.200000, material_new_lambertian(vec3(0.422673, 0.922647, 0.112488))));
// world_add(sphere_new(vec3(2.433369, 0.200000, -2.538340), 0.200000, material_new_lambertian(vec3(0.027875, 0.007925, 0.103471))));
// world_add(sphere_new(vec3(2.476327, 0.200000, -1.222858), 0.200000, material_new_lambertian(vec3(0.298070, 0.155797, 0.696187))));
// world_add(sphere_new(vec3(2.257143, 0.200000, -0.494861), 0.200000, material_new_lambertian(vec3(0.392493, 0.898109, 0.619451))));
// world_add(sphere_new(vec3(2.314603, 0.200000, 0.386785), 0.200000, material_new_lambertian(vec3(0.041665, 0.066173, 0.645386))));
// world_add(sphere_new(vec3(2.588253, 0.200000, 1.451415), 0.200000, material_new_lambertian(vec3(0.173666, 0.281586, 0.322262))));
// world_add(sphere_new(vec3(2.189245, 0.200000, 2.473772), 0.200000, material_new_lambertian(vec3(0.476366, 0.062532, 0.099317))));
// world_add(sphere_new(vec3(2.858086, 0.200000, 4.316443), 0.200000, material_new_lambertian(vec3(0.147476, 0.718038, 0.023072))));
// world_add(sphere_new(vec3(2.517966, 0.200000, 5.445674), 0.200000, material_new_lambertian(vec3(0.439383, 0.450896, 0.003753))));
// world_add(sphere_new(vec3(2.770440, 0.200000, 6.086135), 0.200000, material_new_lambertian(vec3(0.474085, 0.036099, 0.033058))));
// world_add(sphere_new(vec3(2.051445, 0.200000, 7.066277), 0.200000, material_new_lambertian(vec3(0.101432, 0.288337, 0.191759))));
// world_add(sphere_new(vec3(2.550377, 0.200000, 8.763298), 0.200000, material_new_lambertian(vec3(0.309649, 0.422398, 0.817702))));
// world_add(sphere_new(vec3(2.245112, 0.200000, 9.474926), 0.200000, material_new_lambertian(vec3(0.134720, 0.030031, 0.435740))));
// world_add(sphere_new(vec3(2.475668, 0.200000, 10.733717), 0.200000, material_new_lambertian(vec3(0.039877, 0.110592, 0.177490))));
// world_add(sphere_new(vec3(3.023511, 0.200000, -10.584045), 0.200000, material_new_lambertian(vec3(0.126814, 0.024840, 0.027225))));
// world_add(sphere_new(vec3(3.609458, 0.200000, -9.583441), 0.200000, material_new_lambertian(vec3(0.600390, 0.399366, 0.293605))));
// world_add(sphere_new(vec3(3.813617, 0.200000, -8.169930), 0.200000, material_new_lambertian(vec3(0.628822, 0.020515, 0.492754))));
// world_add(sphere_new(vec3(3.170074, 0.200000, -7.738188), 0.200000, material_new_lambertian(vec3(0.130937, 0.517106, 0.023983))));
// world_add(sphere_new(vec3(3.761541, 0.200000, -5.968523), 0.200000, material_new_lambertian(vec3(0.012788, 0.126096, 0.213994))));
// world_add(sphere_new(vec3(3.641896, 0.200000, -4.699268), 0.200000, material_new_lambertian(vec3(0.063366, 0.148035, 0.001709))));
// world_add(sphere_new(vec3(3.246348, 0.200000, -3.929521), 0.200000, material_new_lambertian(vec3(0.263834, 0.221955, 0.406613))));
// world_add(sphere_new(vec3(3.402689, 0.200000, -2.311823), 0.200000, material_new_lambertian(vec3(0.452710, 0.013175, 0.078962))));
// world_add(sphere_new(vec3(3.470998, 0.200000, -1.194705), 0.200000, material_new_lambertian(vec3(0.040214, 0.016429, 0.400268))));
// world_add(sphere_new(vec3(3.791232, 0.200000, -0.894858), 0.200000, material_new_lambertian(vec3(0.237402, 0.537127, 0.255115))));
// world_add(sphere_new(vec3(3.772170, 0.200000, 1.180374), 0.200000, material_new_lambertian(vec3(0.152105, 0.081203, 0.683533))));
// world_add(sphere_new(vec3(3.009229, 0.200000, 2.470806), 0.200000, material_new_lambertian(vec3(0.471125, 0.243826, 0.042404))));
// world_add(sphere_new(vec3(3.550322, 0.200000, 3.098495), 0.200000, material_new_lambertian(vec3(0.154461, 0.136873, 0.182875))));
// world_add(sphere_new(vec3(3.489126, 0.200000, 4.625031), 0.200000, material_new_lambertian(vec3(0.230639, 0.126806, 0.145926))));
// world_add(sphere_new(vec3(3.692682, 0.200000, 6.562352), 0.200000, material_new_lambertian(vec3(0.122111, 0.807302, 0.561923))));
// world_add(sphere_new(vec3(3.246486, 0.200000, 7.251540), 0.200000, material_new_lambertian(vec3(0.123075, 0.180830, 0.522134))));
world_add(sphere_new(vec3(3.253352, 0.200000, 8.893957), 0.200000, material_new_lambertian(vec3(0.203737, 0.009413, 0.608415))));
world_add(sphere_new(vec3(3.848115, 0.200000, 9.889892), 0.200000, material_new_lambertian(vec3(0.125007, 0.634484, 0.066238))));
world_add(sphere_new(vec3(4.854598, 0.200000, -10.928147), 0.200000, material_new_lambertian(vec3(0.009655, 0.165256, 0.028057))));
world_add(sphere_new(vec3(4.653459, 0.200000, -7.249556), 0.200000, material_new_lambertian(vec3(0.388326, 0.548381, 0.926932))));
world_add(sphere_new(vec3(4.016480, 0.200000, -6.822950), 0.200000, material_new_lambertian(vec3(0.072834, 0.264129, 0.568185))));
world_add(sphere_new(vec3(4.307517, 0.200000, -5.189816), 0.200000, material_new_lambertian(vec3(0.135767, 0.093583, 0.141522))));
world_add(sphere_new(vec3(4.858855, 0.200000, -4.458110), 0.200000, material_new_lambertian(vec3(0.208107, 0.245657, 0.243189))));
world_add(sphere_new(vec3(4.359264, 0.200000, -3.537956), 0.200000, material_new_lambertian(vec3(0.045056, 0.059387, 0.239307))));
world_add(sphere_new(vec3(4.387884, 0.200000, -1.650569), 0.200000, material_new_lambertian(vec3(0.370554, 0.000754, 0.230108))));
world_add(sphere_new(vec3(4.786178, 0.200000, -0.611786), 0.200000, material_new_lambertian(vec3(0.148103, 0.130540, 0.011947))));
world_add(sphere_new(vec3(4.694632, 0.200000, 0.631980), 0.200000, material_new_lambertian(vec3(0.157548, 0.199185, 0.293093))));
world_add(sphere_new(vec3(4.188861, 0.200000, 1.646620), 0.200000, material_new_lambertian(vec3(0.039444, 0.033988, 0.040209))));
world_add(sphere_new(vec3(4.603360, 0.200000, 2.662578), 0.200000, material_new_lambertian(vec3(0.492647, 0.000074, 0.029823))));
world_add(sphere_new(vec3(4.136070, 0.200000, 3.298865), 0.200000, material_new_lambertian(vec3(0.179209, 0.544888, 0.276357))));
world_add(sphere_new(vec3(4.793210, 0.200000, 4.416587), 0.200000, material_new_lambertian(vec3(0.007261, 0.170077, 0.172194))));
world_add(sphere_new(vec3(4.549361, 0.200000, 5.684936), 0.200000, material_new_lambertian(vec3(0.719029, 0.071592, 0.041731))));
world_add(sphere_new(vec3(4.689276, 0.200000, 6.645988), 0.200000, material_new_lambertian(vec3(0.169596, 0.127200, 0.003011))));
// world_add(sphere_new(vec3(4.212784, 0.200000, 7.226270), 0.200000, material_new_lambertian(vec3(0.178894, 0.105178, 0.034959))));
// world_add(sphere_new(vec3(4.691391, 0.200000, 8.073803), 0.200000, material_new_lambertian(vec3(0.003636, 0.587201, 0.117632))));
// world_add(sphere_new(vec3(4.521262, 0.200000, 9.413923), 0.200000, material_new_lambertian(vec3(0.355980, 0.332043, 0.656752))));
// world_add(sphere_new(vec3(4.389450, 0.200000, 10.226820), 0.200000, material_new_lambertian(vec3(0.539622, 0.116324, 0.051618))));
// world_add(sphere_new(vec3(5.617121, 0.200000, -10.732832), 0.200000, material_new_lambertian(vec3(0.307448, 0.171498, 0.092560))));
// world_add(sphere_new(vec3(5.331523, 0.200000, -9.616593), 0.200000, material_new_lambertian(vec3(0.219784, 0.539473, 0.175860))));
// world_add(sphere_new(vec3(5.726301, 0.200000, -7.921940), 0.200000, material_new_lambertian(vec3(0.021111, 0.029359, 0.076351))));
// world_add(sphere_new(vec3(5.523707, 0.200000, -4.118403), 0.200000, material_new_lambertian(vec3(0.056847, 0.010092, 0.699475))));
// world_add(sphere_new(vec3(5.046144, 0.200000, -3.775900), 0.200000, material_new_lambertian(vec3(0.008460, 0.079250, 0.011885))));
// world_add(sphere_new(vec3(5.175869, 0.200000, -2.313196), 0.200000, material_new_lambertian(vec3(0.199315, 0.467275, 0.093837))));
// world_add(sphere_new(vec3(5.104044, 0.200000, -1.107169), 0.200000, material_new_lambertian(vec3(0.078432, 0.044430, 0.211840))));
// world_add(sphere_new(vec3(5.492999, 0.200000, -0.915787), 0.200000, material_new_lambertian(vec3(0.209931, 0.059499, 0.520422))));
// world_add(sphere_new(vec3(5.454299, 0.200000, 1.855449), 0.200000, material_new_lambertian(vec3(0.274277, 0.304049, 0.646258))));
// world_add(sphere_new(vec3(5.152660, 0.200000, 2.150874), 0.200000, material_new_lambertian(vec3(0.020384, 0.057413, 0.361403))));
// world_add(sphere_new(vec3(5.616105, 0.200000, 4.101489), 0.200000, material_new_lambertian(vec3(0.120666, 0.043572, 0.460593))));
// world_add(sphere_new(vec3(5.678234, 0.200000, 6.811393), 0.200000, material_new_lambertian(vec3(0.030839, 0.177295, 0.232614))));
world_add(sphere_new(vec3(5.058477, 0.200000, 7.816666), 0.200000, material_new_lambertian(vec3(0.606542, 0.329558, 0.020104))));
world_add(sphere_new(vec3(5.547411, 0.200000, 8.109372), 0.200000, material_new_lambertian(vec3(0.097559, 0.163104, 0.349226))));
world_add(sphere_new(vec3(5.044716, 0.200000, 9.150407), 0.200000, material_new_lambertian(vec3(0.015133, 0.073300, 0.457788))));
world_add(sphere_new(vec3(5.052983, 0.200000, 10.714710), 0.200000, material_new_lambertian(vec3(0.471844, 0.167640, 0.174701))));
world_add(sphere_new(vec3(6.003928, 0.200000, -10.974264), 0.200000, material_new_lambertian(vec3(0.037304, 0.434891, 0.720521))));
world_add(sphere_new(vec3(6.174248, 0.200000, -9.488433), 0.200000, material_new_lambertian(vec3(0.347510, 0.644033, 0.653901))));
world_add(sphere_new(vec3(6.026423, 0.200000, -8.355275), 0.200000, material_new_lambertian(vec3(0.311465, 0.027858, 0.025730))));
world_add(sphere_new(vec3(6.039909, 0.200000, -7.446547), 0.200000, material_new_lambertian(vec3(0.035671, 0.707263, 0.004974))));
// world_add(sphere_new(vec3(6.115854, 0.200000, -6.915650), 0.200000, material_new_lambertian(vec3(0.475661, 0.001889, 0.042066))));
// world_add(sphere_new(vec3(6.008185, 0.200000, -5.826164), 0.200000, material_new_lambertian(vec3(0.005066, 0.019006, 0.166004))));
// world_add(sphere_new(vec3(6.348524, 0.200000, -2.439763), 0.200000, material_new_lambertian(vec3(0.193187, 0.823454, 0.105944))));
// world_add(sphere_new(vec3(6.374755, 0.200000, -1.377441), 0.200000, material_new_lambertian(vec3(0.248535, 0.056330, 0.314721))));
// world_add(sphere_new(vec3(6.147880, 0.200000, -0.192233), 0.200000, material_new_lambertian(vec3(0.191456, 0.077684, 0.267306))));
// world_add(sphere_new(vec3(6.366515, 0.200000, 0.138844), 0.200000, material_new_lambertian(vec3(0.138988, 0.170251, 0.001380))));
// world_add(sphere_new(vec3(6.575838, 0.200000, 1.751460), 0.200000, material_new_lambertian(vec3(0.155068, 0.874774, 0.728850))));
// world_add(sphere_new(vec3(6.657662, 0.200000, 2.371926), 0.200000, material_new_lambertian(vec3(0.308058, 0.589887, 0.138782))));
// world_add(sphere_new(vec3(6.225501, 0.200000, 3.127940), 0.200000, material_new_lambertian(vec3(0.161546, 0.302615, 0.052685))));
// world_add(sphere_new(vec3(6.602097, 0.200000, 4.107642), 0.200000, material_new_lambertian(vec3(0.073070, 0.273733, 0.217113))));
// world_add(sphere_new(vec3(6.409418, 0.200000, 6.778103), 0.200000, material_new_lambertian(vec3(0.424396, 0.021028, 0.000875))));
// world_add(sphere_new(vec3(6.123408, 0.200000, 7.286312), 0.200000, material_new_lambertian(vec3(0.434469, 0.358797, 0.195731))));
// world_add(sphere_new(vec3(6.618192, 0.200000, 8.372832), 0.200000, material_new_lambertian(vec3(0.197941, 0.337186, 0.144319))));
// world_add(sphere_new(vec3(6.844353, 0.200000, 9.056774), 0.200000, material_new_lambertian(vec3(0.142664, 0.901492, 0.520076))));
// world_add(sphere_new(vec3(6.747862, 0.200000, 10.592703), 0.200000, material_new_lambertian(vec3(0.210195, 0.389496, 0.002857))));
// world_add(sphere_new(vec3(7.355089, 0.200000, -10.220331), 0.200000, material_new_lambertian(vec3(0.004409, 0.172869, 0.133061))));
// world_add(sphere_new(vec3(7.809360, 0.200000, -9.857860), 0.200000, material_new_lambertian(vec3(0.278814, 0.324448, 0.041904))));
// world_add(sphere_new(vec3(7.157714, 0.200000, -8.242689), 0.200000, material_new_lambertian(vec3(0.032755, 0.121344, 0.253532))));
// world_add(sphere_new(vec3(7.315372, 0.200000, -7.771477), 0.200000, material_new_lambertian(vec3(0.727132, 0.111195, 0.330344))));
// world_add(sphere_new(vec3(7.754372, 0.200000, -6.469097), 0.200000, material_new_lambertian(vec3(0.396136, 0.656369, 0.285211))));
// world_add(sphere_new(vec3(7.268020, 0.200000, -3.296826), 0.200000, material_new_lambertian(vec3(0.459511, 0.475719, 0.236688))));
// world_add(sphere_new(vec3(7.292932, 0.200000, -2.825449), 0.200000, material_new_lambertian(vec3(0.325913, 0.388341, 0.301442))));
// world_add(sphere_new(vec3(7.561583, 0.200000, -1.185037), 0.200000, material_new_lambertian(vec3(0.309654, 0.230135, 0.162330))));
// world_add(sphere_new(vec3(7.669253, 0.200000, -0.575503), 0.200000, material_new_lambertian(vec3(0.107268, 0.023862, 0.178511))));
// world_add(sphere_new(vec3(7.758519, 0.200000, 2.367806), 0.200000, material_new_lambertian(vec3(0.411538, 0.528854, 0.157908))));
// world_add(sphere_new(vec3(7.722428, 0.200000, 3.599927), 0.200000, material_new_lambertian(vec3(0.076433, 0.195964, 0.223908))));
// world_add(sphere_new(vec3(7.298453, 0.200000, 5.579217), 0.200000, material_new_lambertian(vec3(0.829967, 0.139123, 0.613312))));
// world_add(sphere_new(vec3(7.624015, 0.200000, 6.513297), 0.200000, material_new_lambertian(vec3(0.164426, 0.697133, 0.088994))));
// world_add(sphere_new(vec3(7.694934, 0.200000, 7.339653), 0.200000, material_new_lambertian(vec3(0.476054, 0.046471, 0.445224))));
// world_add(sphere_new(vec3(7.652773, 0.200000, 8.102039), 0.200000, material_new_lambertian(vec3(0.037552, 0.004446, 0.094365))));
// world_add(sphere_new(vec3(7.898791, 0.200000, 9.808921), 0.200000, material_new_lambertian(vec3(0.319317, 0.437308, 0.009118))));
// world_add(sphere_new(vec3(7.058394, 0.200000, 10.158510), 0.200000, material_new_lambertian(vec3(0.390329, 0.404793, 0.301330))));
// world_add(sphere_new(vec3(8.205258, 0.200000, -10.746702), 0.200000, material_new_lambertian(vec3(0.016682, 0.130120, 0.162584))));
// world_add(sphere_new(vec3(8.222507, 0.200000, -8.948857), 0.200000, material_new_lambertian(vec3(0.072035, 0.676721, 0.182858))));
// world_add(sphere_new(vec3(8.582018, 0.200000, -6.861953), 0.200000, material_new_lambertian(vec3(0.822004, 0.440670, 0.068858))));
// world_add(sphere_new(vec3(8.329490, 0.200000, -5.250133), 0.200000, material_new_lambertian(vec3(0.318736, 0.161179, 0.320728))));
// world_add(sphere_new(vec3(8.701993, 0.200000, -3.796857), 0.200000, material_new_lambertian(vec3(0.140003, 0.239332, 0.399558))));
// world_add(sphere_new(vec3(8.809250, 0.200000, -2.778069), 0.200000, material_new_lambertian(vec3(0.181242, 0.117818, 0.130816))));
// world_add(sphere_new(vec3(8.377502, 0.200000, -1.134910), 0.200000, material_new_lambertian(vec3(0.573631, 0.009494, 0.229500))));
// world_add(sphere_new(vec3(8.344954, 0.200000, -0.472997), 0.200000, material_new_lambertian(vec3(0.087736, 0.049359, 0.232017))));
// world_add(sphere_new(vec3(8.012745, 0.200000, 1.158922), 0.200000, material_new_lambertian(vec3(0.636147, 0.180742, 0.311922))));
// world_add(sphere_new(vec3(8.360967, 0.200000, 2.826939), 0.200000, material_new_lambertian(vec3(0.301271, 0.298951, 0.220090))));
// world_add(sphere_new(vec3(8.514203, 0.200000, 3.211933), 0.200000, material_new_lambertian(vec3(0.199712, 0.101657, 0.003248))));
// world_add(sphere_new(vec3(8.437846, 0.200000, 4.172463), 0.200000, material_new_lambertian(vec3(0.025949, 0.117427, 0.624853))));
// world_add(sphere_new(vec3(8.097287, 0.200000, 5.435511), 0.200000, material_new_lambertian(vec3(0.628501, 0.319914, 0.662289))));
// world_add(sphere_new(vec3(8.060234, 0.200000, 6.529804), 0.200000, material_new_lambertian(vec3(0.031272, 0.081360, 0.146070))));
// world_add(sphere_new(vec3(8.636842, 0.200000, 7.715781), 0.200000, material_new_lambertian(vec3(0.493374, 0.516282, 0.734975))));
// world_add(sphere_new(vec3(8.765770, 0.200000, 9.466411), 0.200000, material_new_lambertian(vec3(0.503445, 0.247663, 0.281448))));
// world_add(sphere_new(vec3(8.797714, 0.200000, 10.121320), 0.200000, material_new_lambertian(vec3(0.121923, 0.287376, 0.319318))));
// world_add(sphere_new(vec3(9.736216, 0.200000, -10.196573), 0.200000, material_new_lambertian(vec3(0.720489, 0.006645, 0.009772))));
// world_add(sphere_new(vec3(9.352123, 0.200000, -9.211515), 0.200000, material_new_lambertian(vec3(0.684419, 0.207199, 0.410926))));
// world_add(sphere_new(vec3(9.803674, 0.200000, -8.289877), 0.200000, material_new_lambertian(vec3(0.480957, 0.374283, 0.629879))));
// world_add(sphere_new(vec3(9.690045, 0.200000, -6.256038), 0.200000, material_new_lambertian(vec3(0.252991, 0.616008, 0.110919))));
// world_add(sphere_new(vec3(9.424085, 0.200000, -5.342943), 0.200000, material_new_lambertian(vec3(0.120768, 0.330910, 0.089742))));
// world_add(sphere_new(vec3(9.813013, 0.200000, -4.676965), 0.200000, material_new_lambertian(vec3(0.932773, 0.679638, 0.313875))));
// world_add(sphere_new(vec3(9.215613, 0.200000, -3.475579), 0.200000, material_new_lambertian(vec3(0.371448, 0.088976, 0.384858))));
// world_add(sphere_new(vec3(9.759481, 0.200000, -2.896698), 0.200000, material_new_lambertian(vec3(0.412531, 0.132314, 0.019700))));
// world_add(sphere_new(vec3(9.117832, 0.200000, -1.923066), 0.200000, material_new_lambertian(vec3(0.379933, 0.087664, 0.472595))));
// world_add(sphere_new(vec3(9.492917, 0.200000, -0.547761), 0.200000, material_new_lambertian(vec3(0.182273, 0.452177, 0.693358))));
// world_add(sphere_new(vec3(9.613303, 0.200000, 0.011728), 0.200000, material_new_lambertian(vec3(0.092947, 0.537273, 0.171686))));
// world_add(sphere_new(vec3(9.056114, 0.200000, 1.175594), 0.200000, material_new_lambertian(vec3(0.653357, 0.383002, 0.012880))));
// world_add(sphere_new(vec3(9.268899, 0.200000, 2.223057), 0.200000, material_new_lambertian(vec3(0.096493, 0.424215, 0.146845))));
// world_add(sphere_new(vec3(9.357479, 0.200000, 3.492999), 0.200000, material_new_lambertian(vec3(0.497388, 0.363370, 0.176109))));
// world_add(sphere_new(vec3(9.087701, 0.200000, 4.189163), 0.200000, material_new_lambertian(vec3(0.063743, 0.086766, 0.103873))));
// world_add(sphere_new(vec3(9.380056, 0.200000, 5.605695), 0.200000, material_new_lambertian(vec3(0.293758, 0.139652, 0.055501))));
// world_add(sphere_new(vec3(9.023017, 0.200000, 6.471932), 0.200000, material_new_lambertian(vec3(0.185448, 0.208607, 0.374940))));
// world_add(sphere_new(vec3(9.854295, 0.200000, 7.409748), 0.200000, material_new_lambertian(vec3(0.148198, 0.549944, 0.159839))));
// world_add(sphere_new(vec3(9.862481, 0.200000, 9.574712), 0.200000, material_new_lambertian(vec3(0.765350, 0.008606, 0.235799))));
// world_add(sphere_new(vec3(9.386840, 0.200000, 10.826911), 0.200000, material_new_lambertian(vec3(0.089196, 0.487001, 0.045594))));
// world_add(sphere_new(vec3(10.174056, 0.200000, -10.128263), 0.200000, material_new_lambertian(vec3(0.005711, 0.473478, 0.413344))));
// world_add(sphere_new(vec3(10.872231, 0.200000, -9.900790), 0.200000, material_new_lambertian(vec3(0.305351, 0.335677, 0.067792))));
// world_add(sphere_new(vec3(10.823341, 0.200000, -8.112690), 0.200000, material_new_lambertian(vec3(0.539657, 0.207191, 0.003156))));
// world_add(sphere_new(vec3(10.334901, 0.200000, -7.889831), 0.200000, material_new_lambertian(vec3(0.003655, 0.051635, 0.149474))));
// world_add(sphere_new(vec3(10.418784, 0.200000, -6.461873), 0.200000, material_new_lambertian(vec3(0.288085, 0.073273, 0.500898))));
// world_add(sphere_new(vec3(10.129258, 0.200000, -5.674136), 0.200000, material_new_lambertian(vec3(0.152473, 0.461963, 0.606719))));
// world_add(sphere_new(vec3(10.839161, 0.200000, -4.514609), 0.200000, material_new_lambertian(vec3(0.132318, 0.707132, 0.411638))));
world_add(sphere_new(vec3(10.480667, 0.200000, -3.924604), 0.200000, material_new_lambertian(vec3(0.526081, 0.005696, 0.030672))));
world_add(sphere_new(vec3(10.815787, 0.200000, -1.133811), 0.200000, material_new_lambertian(vec3(0.068619, 0.202986, 0.454007))));
world_add(sphere_new(vec3(10.731822, 0.200000, 0.288977), 0.200000, material_new_lambertian(vec3(0.343254, 0.020722, 0.206464))));
world_add(sphere_new(vec3(10.574355, 0.200000, 1.303204), 0.200000, material_new_lambertian(vec3(0.331582, 0.157878, 0.039533))));
world_add(sphere_new(vec3(10.205231, 0.200000, 3.275436), 0.200000, material_new_lambertian(vec3(0.696670, 0.105404, 0.183944))));
world_add(sphere_new(vec3(10.653240, 0.200000, 5.243959), 0.200000, material_new_lambertian(vec3(0.396251, 0.014514, 0.088905))));
world_add(sphere_new(vec3(10.508545, 0.200000, 6.687765), 0.200000, material_new_lambertian(vec3(0.588092, 0.060319, 0.368795))));
world_add(sphere_new(vec3(10.093771, 0.200000, 7.706113), 0.200000, material_new_lambertian(vec3(0.073221, 0.166891, 0.080368))));
world_add(sphere_new(vec3(10.313422, 0.200000, 9.543538), 0.200000, material_new_lambertian(vec3(0.024323, 0.309233, 0.115325))));
world_add(sphere_new(vec3(10.550514, 0.200000, 10.269008), 0.200000, material_new_lambertian(vec3(0.130278, 0.086473, 0.043980))));
}

void random_scene() {
    world_add(sphere_new(vec3(0.0, -1000.0, 0.0), 1000.0, material_new_lambertian(vec3(0.5, 0.5, 0.5))));

    
    random_scene1();
    random_scene2();
    random_scene3();
    world_add(sphere_new(vec3(0, 1, 0), 1.0, material_new_dielectric(1.5)));

    world_add(
        sphere_new(vec3(-4, 1, 0), 1.0, material_new_lambertian(vec3(0.4, 0.2, 0.1))));

    world_add(
        sphere_new(vec3(4, 1, 0), 1.0, material_new_metal(vec3(0.7, 0.6, 0.5), 0.0)));
}

void simple_scene() {
    world_add(sphere_new(vec3(0.0, -100.5,-1.0), 100.0, material_new_lambertian(vec3(0.8, 0.8, 0.0))));
    world_add(sphere_new(vec3(0.0, 0.0,-1.0), 0.5, material_new_lambertian(vec3(0.1, 0.2, 0.5))));
    world_add(sphere_new(vec3(1.0, 0.0,-1.0), 0.5, material_new_metal(vec3(0.7, 0.6, 0.5), 0.0)));
    world_add(sphere_new(vec3(-1.0, 0.0,-1.0), 0.5, material_new_dielectric(1.5)));
    world_add(sphere_new(vec3(-1.0, 0.0,-1.0), -0.45, material_new_dielectric(1.5)));
}

void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

    // simple_scene();
    random_scene();

    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    vec4 prev = texture(iChannel0, uv);
    vec3 prevLinear = toLinear(prev.xyz);
    prevLinear *= prev.w;

    uv = (gl_FragCoord.xy + hash2(gSeed)) / iResolution.xy;
    vec3 lookfrom = vec3(13.0,2.0,3.0);
    vec3 lookat = vec3(0.0,0.0,0.0);
    float aspect_ratio = iResolution.x / iResolution.y;
    float dist_to_focus = 10.0;//length((lookfrom - lookat));
    camera c = camera_new(lookfrom, lookat, vec3(0.0, 1.0, 0.0), 20.0, aspect_ratio, 0.1, dist_to_focus);
    ray r = camera_get_ray(c, uv);
    vec3 color = ray_color_depth(r, 50);

    if (prev.w > 500.0) {
        gl_FragColor = prev;
        return;
    }

    color = (color + prevLinear);
    float w = prev.w + 1.0;
    color = color / w;
    color = toGamma(color);

    gl_FragColor = vec4(color, w);
}

