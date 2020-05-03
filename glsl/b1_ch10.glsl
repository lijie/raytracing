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
};

camera camera_new(vec3 lookfrom, vec3 lookat, vec3 vup, float fov, float aspect) {
    camera c;

    c.origin = lookfrom;

    float half_height = tan(radians(fov) / 2.0);
    float half_width = aspect * half_height;

    vec3 w = normalize(lookfrom - lookat);
    vec3 u = normalize(cross(vup, w));
    vec3 v = cross(w, u);

    c.lower_left_corner = c.origin - half_width * u - half_height * v - w;
    c.horizontal = 2.0 * half_width * u;
    c.vertical = 2.0 * half_height * v;
    return c;
}

ray camera_get_ray(camera c, vec2 uv) {
    vec3 direction = normalize(c.lower_left_corner + uv.x * c.horizontal + uv.y * c.vertical - c.origin);
    return ray_new(c.origin, direction);
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

float gSeed = 1.0;
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

void main()
{
    gSeed = float(baseHash(floatBitsToUint(gl_FragCoord.xy))) / float(0xffffffffU) + iTime;

    world_add(sphere_new(vec3(0.0, -100.5,-1.0), 100.0, material_new_lambertian(vec3(0.8, 0.8, 0.0))));
    world_add(sphere_new(vec3(0.0, 0.0,-1.0), 0.5, material_new_lambertian(vec3(0.1, 0.2, 0.5))));
    world_add(sphere_new(vec3(1.0, 0.0,-1.0), 0.5, material_new_metal(vec3(0.7, 0.6, 0.5), 0.0)));
    world_add(sphere_new(vec3(-1.0, 0.0,-1.0), 0.5, material_new_dielectric(1.5)));
    world_add(sphere_new(vec3(-1.0, 0.0,-1.0), -0.45, material_new_dielectric(1.5)));

    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    vec4 prev = texture(iChannel0, uv);
    vec3 prevLinear = toLinear(prev.xyz);
    prevLinear *= prev.w;

    uv = (gl_FragCoord.xy + hash2(gSeed)) / iResolution.xy;
    float aspect_ratio = iResolution.x / iResolution.y;
    camera c = camera_new(vec3(-2.0,2.0,1.0), vec3(0.0,0.0,-1.0), vec3(0.0, 1.0, 0.0), 45.0, aspect_ratio);
    ray r = camera_get_ray(c, uv);
    vec3 color = ray_color_depth(r, 50);

    if(prev.w > 500.0)
    {
        gl_FragColor = prev;
        return;
    }

    color = (color + prevLinear);
    float w = prev.w + 1.0;
    color = color / w;
    color = toGamma(color);

    gl_FragColor = vec4(color, w);
}

