#ifndef MC_VOXEL_CAMERA_HPP__
#define MC_VOXEL_CAMERA_HPP__

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "ray.hpp"

namespace mcvoxel
{

/// @brief A camera wraps up a 4x4 matrix which transforms between camera-centred co-ordinates and world co-ordinates.
///
/// Camera co-ordinates have the optical centre of the camera at (0,0,0) and the centre of the plane at (0,0,-f).
/// Positive x is toward the right of the image plane and positive y is toward the top.
///
/// The camera is parameterised in terms of three vectors: the camera centre, the 'up' direction (+ve camera y-axis) and
/// the 'look at' direction (-ve camera z-axis). Moving along the 'look at' axis moves <em>into</em> the scene. Moving
/// along the 'up' direction moves, subjectively for the camera, upwards. The 'up' direction (and derived 'right'
/// direction) have magnitude of one pixel. The 'look at', 'up' and 'right' vectors all have unit magnitude.
///
/// There is a derived direction, 'right', which is the cross-product of 'look at' and 'up'.
///
/// Internally the camera is represented as a transform <em>to</em> world co-ordinates <em>from</em> camera-centric
/// co-ordinates. I.e., the camera centre is translated to the origin and then the co-ordinate system is rotated such
/// that 'look at' lies along the -ve z-axis and 'up' lies along the +ve y-axis.
class camera : protected Eigen::Transform<float, 3, Eigen::AffineCompact>
{
public:
	typedef Eigen::Transform<float, 3, Eigen::AffineCompact> transform_type;

	/// @brief Get a reference to this camera as the underlying Eigen::Transform instance.
	///
	/// This is the transformation <em>from</em> camera-space <em>to</em> world-space.
	const transform_type& transform() const { return *this; }

	/// @brief Return the inverse transform for this camera.
	///
	/// This is the transform which maps <em>from</em> world-space <em>to</em> camera-space.
	const transform_type inverse() const { return transform().inverse(); }

	/// @brief Write the four vectors describing the frame for the camera into the four references provided.
	///
	/// @param centre
	/// @param look_at
	/// @param up
	/// @param right
	void get_frame(Eigen::Vector3f& centre,
		       Eigen::Vector3f& look_at,
		       Eigen::Vector3f& up,
		       Eigen::Vector3f& right) const;

	/// @brief Get the centre of the camera in world-space.
	///
	/// @note It is more efficient to use get_frame() if you want to get multiple frame vectors.
	Eigen::Vector3f centre() const { return transform() * Eigen::Vector3f(0.f,0.f,0.f); }

	/// @brief Set the camera centre in world-space.
	///
	/// @param centre
	void set_centre(Eigen::Vector3f centre);

	/// @brief Get the 'look at' vector of the camera in world-space.
	///
	/// @note It is more efficient to use get_frame() if you want to get multiple frame vectors.
	Eigen::Vector3f look_at() const { return transform() * Eigen::Vector3f(0.f,0.f,-1.f) - centre(); }

	/// @brief Get the 'up' vector of the camera in world-space.
	///
	/// @note It is more efficient to use get_frame() if you want to get multiple frame vectors.
	Eigen::Vector3f up() const { return transform() * Eigen::Vector3f(0.f,1.f,0.f) - centre(); }

	/// @brief Get the 'right' vector of the camera in world-space.
	///
	/// @note It is more efficient to use get_frame() if you want to get multiple frame vectors.
	Eigen::Vector3f right() const { return transform() * Eigen::Vector3f(1.f,0.f,0.f) - centre(); }

	/// @brief Get the focal length of the camera. This is simply the magnitude of the 'look at' vector.
	float focal_length() const { return focal_length_; }

	/// @brief Set the focal length.
	///
	/// This is the magnitude of the 'look at' vector.
	///
	/// @param focal_length
	void set_focal_length(float focal_length) { focal_length_ = focal_length; }

	/// @brief Apply a rotation to the camera.
	///
	/// This is a subjective rotation of the camera. For example, if given a -ve rotation about the +ve x-axis the
	/// camera will subjectively pan <em>down</em>. Similarly a +ve rotation about the -ve z-axis will cause the
	/// camera to subjectively roll clockwise causing the <em>image</em> to roll counter-clockwise.
	///
	/// @tparam RotationType
	/// @param r
	template<typename RotationType>
	void rotate(const RotationType& r) { transform_type::rotate(r); }

	/// @brief Translate the camera co-ordinate system.
	///
	/// This is a subjective translation of the co-ordinate system. This has the effect of translating the camera
	/// centre in world-space by -ve \p t. Note the sign!
	///
	/// @tparam OtherDerived
	/// @param t
	template<typename OtherDerived>
	void translate(const Eigen::MatrixBase<OtherDerived> & t) { transform_type::pretranslate(t); }

	/// @brief Roll the camera clockwise around the camera axis.
	///
	/// Rolling the camera clockwise will cause the <em>image</em> to rotate counter-clockwise.
	///
	/// @param angle
	void roll_clockwise(float angle)
       	{ rotate(Eigen::AngleAxis<float>(angle, Eigen::Vector3f(0.f,0.f,-1.f))); }

	/// @brief Pitch the camera up around the right axis.
	///
	/// Pitching the camera in a +ve direction will cause the <em>image</em> to move down.
	///
	/// @param angle
	void pitch_up(float angle)
       	{ rotate(Eigen::AngleAxis<float>(angle, Eigen::Vector3f(1.f,0.f,0.f))); }

	/// @brief Yaw the camera left around the up acis.
	///
	/// Yawing the camera in a +ve direction will cause the <em>image</em> to move to the right.
	///
	/// @param angle
	void yaw_left(float angle)
       	{ rotate(Eigen::AngleAxis<float>(angle, Eigen::Vector3f(0.f,1.f,0.f))); }

	/// @brief Return an eye-ray for the pixel (\p x, \p y).
	///
	/// @param x
	/// @param y
	///
	/// @return 
	ray eye_ray(float x, float y) const { return ray(centre(), focal_length() * look_at() + x * right() + y * up()); }

	/// @brief Create a camera looking in the (0,0,-1) direction with an up direction of (0,1,0) and unit focal
	/// length.
	camera()
	       	: transform_type(transform_type::Identity())
		, focal_length_(1.f)
	{ }

	camera(const camera& c)
	       	: transform_type(c)
		, focal_length_(c.focal_length_)
	{ }

	const camera& operator = (const camera& c)
	{
		transform_type::operator = (c);
		focal_length_ = c.focal_length_;
		return *this;
	}

protected:
	float focal_length_;
};

}

#endif // MC_VOXEL_CAMERA_HPP__
